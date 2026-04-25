import pandas as pd
import requests
import urllib.parse
import time
import json
import os
import sys
from functools import lru_cache
from datetime import datetime

# =============================================================================
# 🧬 全RNA靶点药物重定位自动化分析引擎
# 支持两种模式：全RNA分析 / 抗菌专项分析
# 完全兼容GitHub Actions v5工作流
# =============================================================================

# --------------------------
# 配置与缓存系统
# --------------------------
SMILES_CACHE_FILE = "smiles_cache.json"
DRUG_SEARCH_CACHE_FILE = "drug_search_cache.json"
DRUG_INDICATION_CACHE_FILE = "drug_indication_cache.json"  # 新增适应症缓存
DRUG_TARGET_CACHE_FILE = "drug_target_cache.json"  # 新增靶点缓存

# 输出目录（与GitHub Actions工作流完全一致）
ALL_RNA_OUTPUT_DIR = "all_rna_drug_discovery_results"
ANTIBACTERIAL_OUTPUT_DIR = "antibacterial_results"

# ... 保留原有ATC/适应症映射字典（ATC_FULL_MAP/ATC_CLASS_MAP/INDICATION_CN_MAP/MANUAL_DRUG_INDICATION） ...

# --------------------------
# 新增：药物适应症补全 & 官方靶点获取函数
# --------------------------
@lru_cache(maxsize=500)
def get_drug_indication_by_name(drug_name):
    """通过药物名称补全适应症（优先ChEMBL，其次PubChem）"""
    if not drug_name:
        return "暂无明确适应症"
    
    drug_name_clean = drug_name.strip().upper()
    cache = load_cache(DRUG_INDICATION_CACHE_FILE)
    
    # 1. 先查缓存
    if drug_name_clean in cache:
        return cache[drug_name_clean]
    
    # 2. 查手动映射表
    if drug_name_clean in MANUAL_DRUG_INDICATION:
        cache[drug_name_clean] = MANUAL_DRUG_INDICATION[drug_name_clean]
        save_cache(cache, DRUG_INDICATION_CACHE_FILE)
        return MANUAL_DRUG_INDICATION[drug_name_clean]
    
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
        "Accept": "application/json"
    }
    indication = "暂无明确适应症"
    max_retries = 2
    retry_delay = 1

    # 3. 调用ChEMBL API按药物名称查适应症
    try:
        safe_name = urllib.parse.quote(drug_name_clean)
        for attempt in range(max_retries):
            try:
                # 先查药物ID
                url = f"https://www.ebi.ac.uk/chembl/api/data/drug.json?pref_name={safe_name}&limit=1"
                r = requests.get(url, headers=headers, timeout=10)
                if r.status_code == 200 and r.json().get('drugs'):
                    drug_id = r.json()['drugs'][0]['drug_chembl_id']
                    # 查药物适应症
                    url2 = f"https://www.ebi.ac.uk/chembl/api/data/drug_indication.json?drug_chembl_id={drug_id}&limit=10"
                    r2 = requests.get(url2, headers=headers, timeout=10)
                    if r2.status_code == 200 and r2.json().get('drug_indications'):
                        indications = []
                        for ind in r2.json()['drug_indications']:
                            if ind.get('indication_description'):
                                indications.append(ind['indication_description'])
                        if indications:
                            indication = "、".join(list(set(indications)))
                            break
            except Exception as e:
                time.sleep(retry_delay)
                continue
    except:
        pass

    # 4. ChEMBL查不到则查PubChem
    if indication == "暂无明确适应症":
        try:
            safe_name = urllib.parse.quote(drug_name_clean)
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{safe_name}/description/JSON"
            r = requests.get(url, headers=headers, timeout=10)
            if r.status_code == 200 and r.json().get('InformationList', {}).get('Information'):
                desc = r.json()['InformationList']['Information'][0].get('Description', '')
                if desc and "indication" in desc.lower():
                    # 提取适应症相关内容
                    import re
                    ind_match = re.search(r'indication[s]?: (.+?)(\.|;|$)', desc, re.IGNORECASE)
                    if ind_match:
                        indication = ind_match.group(1).strip()
        except:
            pass

    # 5. 缓存结果
    cache[drug_name_clean] = indication if indication != "" else "暂无明确适应症"
    save_cache(cache, DRUG_INDICATION_CACHE_FILE)
    
    return cache[drug_name_clean]

@lru_cache(maxsize=500)
def get_drug_targets_by_name(drug_name):
    """获取上市药物的官方靶点（ChEMBL API）"""
    if not drug_name:
        return "暂无靶点信息"
    
    drug_name_clean = drug_name.strip().upper()
    cache = load_cache(DRUG_TARGET_CACHE_FILE)
    
    # 1. 先查缓存
    if drug_name_clean in cache:
        return cache[drug_name_clean]
    
    headers = {
        "User-Agent": "Mozilla/5.0",
        "Accept": "application/json"
    }
    targets = "暂无靶点信息"
    max_retries = 2
    retry_delay = 1

    # 2. 调用ChEMBL API查靶点
    try:
        safe_name = urllib.parse.quote(drug_name_clean)
        for attempt in range(max_retries):
            try:
                # 先查药物ID
                url = f"https://www.ebi.ac.uk/chembl/api/data/drug.json?pref_name={safe_name}&limit=1"
                r = requests.get(url, headers=headers, timeout=10)
                if r.status_code == 200 and r.json().get('drugs'):
                    drug_id = r.json()['drugs'][0]['drug_chembl_id']
                    # 查药物-靶点关联
                    url2 = f"https://www.ebi.ac.uk/chembl/api/data/drug_target.json?drug_chembl_id={drug_id}&limit=10"
                    r2 = requests.get(url2, headers=headers, timeout=10)
                    if r2.status_code == 200 and r2.json().get('drug_targets'):
                        target_list = []
                        for dt in r2.json()['drug_targets']:
                            if dt.get('target_chembl_id'):
                                # 查靶点详情
                                target_id = dt['target_chembl_id']
                                url3 = f"https://www.ebi.ac.uk/chembl/api/data/target.json?target_chembl_id={target_id}&limit=1"
                                r3 = requests.get(url3, headers=headers, timeout=10)
                                if r3.status_code == 200 and r3.json().get('targets'):
                                    target_name = r3.json()['targets'][0].get('pref_name', '')
                                    if target_name:
                                        target_list.append(target_name)
                        if target_list:
                            targets = "、".join(list(set(target_list)))
                            break
            except Exception as e:
                time.sleep(retry_delay)
                continue
    except:
        pass

    # 3. 缓存结果
    cache[drug_name_clean] = targets if targets != "" else "暂无靶点信息"
    save_cache(cache, DRUG_TARGET_CACHE_FILE)
    
    return cache[drug_name_clean]

# --------------------------
# 保留原有通用工具函数（load_cache/save_cache）...
# --------------------------
def load_cache(cache_file):
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except:
            return {}
    return {}

def save_cache(cache_dict, cache_file):
    try:
        with open(cache_file, 'w', encoding='utf-8') as f:
            json.dump(cache_dict, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"⚠️ 缓存保存失败: {e}")

# --------------------------
# 保留原有核心API函数（get_smiles_by_id）...
# --------------------------
@lru_cache(maxsize=1000)
def get_smiles_by_id(ligand_id):
    # 原有逻辑保持不变...
    if not ligand_id or ligand_id == "ZZZ":
        return None
    
    ligand_id = ligand_id.strip().upper()
    cache = load_cache(SMILES_CACHE_FILE)
    if ligand_id in cache:
        return cache[ligand_id]

    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
        "Accept": "application/json"
    }
    
    max_retries = 2
    retry_delay = 1
    smiles_result = None

    # 1. PDBe v2 API
    for attempt in range(max_retries):
        try:
            url = f"https://www.ebi.ac.uk/pdbe/api/v2/compound/summary/{ligand_id}"
            r = requests.get(url, headers=headers, timeout=10)
            if r.status_code == 200:
                data = r.json()
                if ligand_id in data:
                    comp = data[ligand_id][0]
                    if "smiles" in comp:
                        for s in comp["smiles"]:
                            if s.get("name") == "canonical":
                                smiles_result = s.get("value")
                                break
                        if not smiles_result and comp["smiles"]:
                            smiles_result = comp["smiles"][0].get("value")
                        if smiles_result: break
            break
        except: time.sleep(retry_delay)

    # 2. RCSB PDB API
    if not smiles_result:
        for attempt in range(max_retries):
            try:
                url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
                r = requests.get(url, headers=headers, timeout=10)
                if r.status_code == 200:
                    data = r.json()
                    if "rcsb_chem_comp_descriptor" in data:
                        desc = data["rcsb_chem_comp_descriptor"]
                        if "SMILES_stereo" in desc: smiles_result = desc["SMILES_stereo"]
                        elif "SMILES" in desc: smiles_result = desc["SMILES"]
                        if smiles_result: break
            except: time.sleep(retry_delay)

    # 3. PubChem XRef API
    if not smiles_result:
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/PDB/{ligand_id}/cids/JSON"
            r = requests.get(url, headers=headers, timeout=10)
            if r.status_code == 200:
                cid = r.json()["IdentifierList"]["CID"][0]
                url2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
                r2 = requests.get(url2, headers=headers, timeout=10)
                if r2.status_code == 200:
                    smiles_result = r2.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        except: pass

    if smiles_result:
        cache[ligand_id] = smiles_result
        save_cache(cache, SMILES_CACHE_FILE)
        return smiles_result
    
    return None

# --------------------------
# 修改search_chembl_drugs函数：补全适应症 + 添加靶点列
# --------------------------
@lru_cache(maxsize=500)
def search_chembl_drugs(smiles, similarity_threshold):
    if not smiles: return [], "SMILES为空", 0
    
    cache_key = f"{smiles}_{similarity_threshold}"
    drug_cache = load_cache(DRUG_SEARCH_CACHE_FILE)
    if cache_key in drug_cache:
        return drug_cache[cache_key]["drugs"], drug_cache[cache_key]["msg"], drug_cache[cache_key]["total"]
    
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json?limit=1000&expand=drug_indications"
        r = requests.get(url, headers=headers, timeout=30)
        
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            drugs = []
            
            for m in mols:
                if not (m.get('max_phase') and float(m.get('max_phase')) >= 4.0 and m.get('pref_name') and not m.get('withdrawn_flag')):
                    continue
                
                drug_name = m.get('pref_name', '').upper()
                
                # 1. 智能适应症获取（新增：调用补全函数）
                indication_list = []
                # 先查手动映射
                if drug_name in MANUAL_DRUG_INDICATION:
                    indication_list.append(MANUAL_DRUG_INDICATION[drug_name])
                # 再查ChEMBL/PubChem API补全
                else:
                    api_indication = get_drug_indication_by_name(drug_name)
                    if api_indication != "暂无明确适应症":
                        indication_list.append(api_indication)
                # 最后查原有indication映射
                if not indication_list:
                    drug_indications = m.get('drug_indications', [])
                    for ind in drug_indications:
                        mesh_term = ind.get('mesh_heading', '').lower()
                        efo_term = ind.get('efo_term', '').lower()
                        for eng, cn in INDICATION_CN_MAP.items():
                            if eng in mesh_term or eng in efo_term:
                                indication_list.append(cn)
                                break
                # 兜底：再次调用API补全
                if not indication_list:
                    base_name = drug_name.replace(' SULFATE', '').replace(' HYDROCHLORIDE', '')
                    api_indication = get_drug_indication_by_name(base_name)
                    indication_list.append(api_indication if api_indication else "暂无明确适应症")
                
                # 去重并格式化适应症
                final_indication = "、".join(list(set(indication_list)))
                
                # 2. 药物分类
                drug_class = "未分类"
                atc_list = m.get('atc_classifications', [])
                if atc_list:
                    atc_code = atc_list[0][:3]
                    if atc_code in ATC_FULL_MAP:
                        drug_class = ATC_FULL_MAP[atc_code]
                    else:
                        atc_first = atc_list[0][0].upper()
                        if atc_first in ATC_CLASS_MAP:
                            drug_class = ATC_CLASS_MAP[atc_first]
                
                # 3. 获取官方靶点（新增）
                official_targets = get_drug_targets_by_name(drug_name)
                
                drugs.append({
                    "药物名称": m.get('pref_name'),
                    "ChEMBL ID": m.get('molecule_chembl_id'),
                    "相似度(%)": round(float(m.get('similarity', 0)), 2),
                    "分子量": round(float(m.get('molecule_properties', {}).get('full_mwt', 0)), 2),
                    "药物分类": drug_class,
                    "治疗适应症": final_indication,
                    "上市药物官方靶点": official_targets  # 新增列
                })
            
            drugs.sort(key=lambda x: x["相似度(%)"], reverse=True)
            
            drug_cache[cache_key] = {"drugs": drugs, "msg": "Success", "total": len(mols)}
            save_cache(drug_cache, DRUG_SEARCH_CACHE_FILE)
            
            return drugs, "Success", len(mols)
        else:
            return [], f"接口失败 ({r.status_code})", 0
    except Exception as e:
        return [], f"异常: {str(e)}", 0

# --------------------------
# 保留原有靶点加载函数（load_all_rna_targets/load_antibacterial_targets）...
# --------------------------
def load_all_rna_targets():
    """加载所有RNA靶点"""
    print("🔍 正在加载所有RNA靶点...")
    
    try:
        df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    except FileNotFoundError:
        print("❌ 错误：未找到数据文件 'PDB_Dataset_Info_Full.xlsx'")
        exit(1)
    
    # RNA分类
    def categorize(desc):
        desc_lower = str(desc).lower()
        if 'riboswitch' in desc_lower: return "核糖开关 (Riboswitch)"
        elif 'aptamer' in desc_lower: return "适配体 (Aptamer)"
        elif any(word in desc_lower for word in ['quadruplex', 'g-4', 'g4']): return "G-四联体 (G-quadruplex)"
        elif any(word in desc_lower for word in ['ribosomal', 'ribosome', 'rrna']): return "核糖体 (rRNA)"
        elif 'ribozyme' in desc_lower: return "核酶 (Ribozyme)"
        elif any(word in desc_lower for word in ['ires', 'hairpin', 'stem-loop', 'pseudoknot']): return "特殊结构基元 (Special/Motifs)"
        else: return "其他 RNA (Others)"
    
    df['Category'] = df['Description (描述)'].apply(categorize)
    
    # 提取核心配体ID
    def get_sort_id(ligand_text):
        if ligand_text == "No ligands" or str(ligand_text) == "nan": return "ZZZ"
        all_l = str(ligand_text).split(' | ')
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    
    # 筛选有配体的靶点
    valid_df = df[df['MainLigandID'] != "ZZZ"].copy()
    
    print(f"✅ 加载完成：共 {len(df)} 个PDB结构")
    print(f"   - 有核心配体的靶点: {len(valid_df)} 个")
    print(f"   - 无配体的靶点: {len(df)-len(valid_df)} 个（已跳过）")
    
    # 按类别统计
    category_stats = valid_df['Category'].value_counts()
    print(f"\n📊 靶点类别分布:")
    for cat, count in category_stats.items():
        print(f"   - {cat}: {count} 个")
    
    return valid_df

def load_antibacterial_targets():
    """加载抗菌RNA靶点（核糖开关和核糖体）"""
    print("🔍 正在加载抗菌RNA靶点...")
    
    try:
        df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    except FileNotFoundError:
        print("❌ 错误：未找到数据文件 'PDB_Dataset_Info_Full.xlsx'")
        exit(1)
    
    # RNA分类
    def categorize(desc):
        desc_lower = str(desc).lower()
        if 'riboswitch' in desc_lower: return "核糖开关 (Riboswitch)"
        elif any(word in desc_lower for word in ['ribosomal', 'ribosome', 'rrna']): return "核糖体 (rRNA)"
        else: return "其他"
    
    df['Category'] = df['Description (描述)'].apply(categorize)
    
    # 提取核心配体ID
    def get_sort_id(ligand_text):
        if ligand_text == "No ligands" or str(ligand_text) == "nan": return "ZZZ"
        all_l = str(ligand_text).split(' | ')
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    
    # 筛选抗菌靶点
    antibacterial_df = df[
        (df['Category'].isin(["核糖开关 (Riboswitch)", "核糖体 (rRNA)"])) & 
        (df['MainLigandID'] != "ZZZ")
    ].copy()
    
    # 优先筛选细菌相关靶点
    def is_bacterial(desc):
        desc_lower = str(desc).lower()
        return any(word in desc_lower for word in ['bacterial', 'bacteria', 'escherichia', 'coli', 'staphylococcus', 'aureus', 'streptococcus', 'mycobacterium', 'tuberculosis'])
    
    antibacterial_df['IsBacterial'] = antibacterial_df['Description (描述)'].apply(is_bacterial)
    antibacterial_df = antibacterial_df.sort_values(by=['IsBacterial', 'Category'], ascending=[False, True])
    
    print(f"✅ 筛选完成：共找到 {len(antibacterial_df)} 个抗菌RNA靶点")
    print(f"   - 核糖开关: {len(antibacterial_df[antibacterial_df['Category'] == '核糖开关 (Riboswitch)'])} 个")
    print(f"   - 核糖体: {len(antibacterial_df[antibacterial_df['Category'] == '核糖体 (rRNA)'])} 个")
    print(f"   - 明确细菌相关: {len(antibacterial_df[antibacterial_df['IsBacterial']])} 个")
    
    return antibacterial_df

# --------------------------
# 修改batch_analyze_targets函数：添加「上市药物官方靶点」列
# --------------------------
def batch_analyze_targets(targets_df, similarity_threshold, output_dir, report_title, report_subtitle):
    """批量分析靶点并生成报告"""
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n🚀 开始批量药物重定位分析（相似度阈值: {similarity_threshold}%）...")
    
    all_results = []
    total_targets = len(targets_df)
    
    for idx, (_, row) in enumerate(targets_df.iterrows()):
        pdb_id = row['PDB ID']
        ligand_id = row['MainLigandID']
        category = row['Category']
        description = row['Description (描述)']
        
        # 每10个靶点打印一次进度
        if (idx + 1) % 10 == 0 or idx == 0:
            print(f"   正在处理 [{idx+1}/{total_targets}] {pdb_id} (配体: {ligand_id})")
        
        smiles = get_smiles_by_id(ligand_id)
        if not smiles:
            continue
        
        drugs, msg, total = search_chembl_drugs(smiles, similarity_threshold)
        
        for drug in drugs:
            drug["RNA类别"] = category
            drug["PDB ID"] = pdb_id
            drug["配体ID"] = ligand_id
            drug["RNA结构描述"] = description
            all_results.append(drug)
    
    print(f"\n✅ 批量分析完成！共匹配到 {len(all_results)} 条药物记录")
    
    if len(all_results) == 0:
        print("❌ 未匹配到任何药物")
        return
    
    # 优先级排序
    print("\n📊 正在进行药物优先级排序...")
    result_df = pd.DataFrame(all_results)
    
    # 计算优先级分数
    def calculate_priority(row):
        score = 0
        # 相似度权重：60%
        score += row["相似度(%)"] * 0.6
        # 药物分类权重：20%
        priority_classes = ["全身用抗菌药", "抗肿瘤药", "抗病毒药", "心血管系统药物"]
        if row["药物分类"] in priority_classes:
            score += 20
        # 特殊RNA类别权重：20%
        special_categories = ["核糖开关 (Riboswitch)", "核糖体 (rRNA)", "G-四联体 (G-quadruplex)"]
        if row["RNA类别"] in special_categories:
            score += 20
        return round(score, 2)
    
    result_df["优先级分数"] = result_df.apply(calculate_priority, axis=1)
    result_df = result_df.sort_values(by=["优先级分数", "相似度(%)"], ascending=False)
    
    # 重新排列列顺序（新增「上市药物官方靶点」列）
    column_order = [
        "优先级分数", "RNA类别", "PDB ID", "配体ID", "药物名称", 
        "药物分类", "治疗适应症", "上市药物官方靶点",  # 新增列
        "相似度(%)", "分子量", "ChEMBL ID", "RNA结构描述"
    ]
    result_df = result_df[column_order]
    
    # 生成HTML报告
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # 统计信息
    total_targets = len(targets_df)
    total_drugs = len(result_df)
    unique_drugs = result_df["药物名称"].nunique()
    top_20_drugs = result_df.head(20).to_html(index=False, classes="table table-striped table-hover")
    
    # 按RNA类别统计
    category_stats = result_df['RNA类别'].value_counts().to_dict()
    category_stats_html = "".join([f"<li><strong>{cat}</strong>: {count} 条药物记录</li>" for cat, count in category_stats.items()])
    
    # 按药物分类统计
    drug_class_stats = result_df['药物分类'].value_counts().head(10).to_dict()
    drug_class_stats_html = "".join([f"<li><strong>{cls}</strong>: {count} 条</li>" for cls, count in drug_class_stats.items()])
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="zh-CN">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{report_title}</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; }}
            .header {{ background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); color: white; padding: 2rem 0; margin-bottom: 2rem; }}
            .section {{ margin-bottom: 2rem; padding: 1.5rem; background: #f8f9fa; border-radius: 10px; }}
            .stat-card {{ background: white; padding: 1.5rem; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); text-align: center; }}
            .stat-number {{ font-size: 2rem; font-weight: bold; color: #11998e; }}
            table {{ font-size: 0.85rem; }}
            .footer {{ text-align: center; padding: 2rem; color: #6c757d; border-top: 1px solid #eee; margin-top: 3rem; }}
        </style>
    </head>
    <body>
        <div class="header">
            <div class="container">
                <h1 class="text-center">{report_title}</h1>
                <p class="text-center lead">{report_subtitle}</p>
                <p class="text-center">报告生成时间: {timestamp} | 相似度阈值: {similarity_threshold}%</p>
            </div>
        </div>

        <div class="container">
            <div class="section">
                <h2>📊 分析摘要</h2>
                <div class="row">
                    <div class="col-md-3">
                        <div class="stat-card">
                            <div class="stat-number">{total_targets}</div>
                            <div>分析靶点总数</div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="stat-card">
                            <div class="stat-number">{total_drugs}</div>
                            <div>匹配药物总数</div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="stat-card">
                            <div class="stat-number">{unique_drugs}</div>
                            <div>唯一药物数量</div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="stat-card">
                            <div class="stat-number">{len(category_stats)}</div>
                            <div>RNA类别数量</div>
                        </div>
                    </div>
                </div>
            </div>

            <div class="section">
                <h2>📈 数据分布统计</h2>
                <div class="row">
                    <div class="col-md-6">
                        <h5>RNA类别药物分布</h5>
                        <ul>{category_stats_html}</ul>
                    </div>
                    <div class="col-md-6">
                        <h5>药物分类TOP 10</h5>
                        <ul>{drug_class_stats_html}</ul>
                    </div>
                </div>
            </div>

            <div class="section">
                <h2>🏆 高优先级候选药物 TOP 20</h2>
                <p>基于相似度、药物分类和RNA类别综合排序</p>
                {top_20_drugs}
            </div>

            <div class="section">
                <h2>💡 研究建议</h2>
                <div class="row">
                    <div class="col-md-6">
                        <h5>🔬 抗菌药物方向</h5>
                        <ul>
                            <li>优先关注核糖开关和核糖体靶点</li>
                            <li>测试氨基糖苷类、大环内酯类药物的抗菌活性</li>
                            <li>针对超级细菌开发新型RNA靶向药物</li>
                        </ul>
                    </div>
                    <div class="col-md-6">
                        <h5>🎯 抗肿瘤药物方向</h5>
                        <ul>
                            <li>重点研究G-四联体靶点</li>
                            <li>筛选能够稳定G-四联体结构的药物</li>
                            <li>探索老药新用于肿瘤治疗的可能性</li>
                        </ul>
                    </div>
                </div>
            </div>
        </div>

        <div class="footer">
            <p>本报告由RNA结构精细化分类与AIDD药物重定位平台自动生成</p>
            <p>© 2025 RNA药物重定位自动化系统</p>
        </div>
    </body>
    </html>
    """
    
    report_path = os.path.join(output_dir, "drug_discovery_report.html")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    csv_path = os.path.join(output_dir, "candidate_drugs.csv")
    result_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    
    print(f"\n📄 报告生成完成！")
    print(f"   - 完整报告: {report_path}")
    print(f"   - 候选药物列表: {csv_path}")
    
    # 打印TOP 5
    print(f"\n📌 高优先级候选药物TOP 5:")
    for i in range(min(5, len(result_df))):
        row = result_df.iloc[i]
        print(f"   {i+1}. {row['药物名称']} (相似度: {row['相似度(%)']}%, RNA类别: {row['RNA类别']}, 靶点: {row['PDB ID']}, 官方靶点: {row['上市药物官方靶点']})")
    
    return result_df

# --------------------------
# 主函数保持不变（已支持传入阈值参数）
# --------------------------
def main():
    # 解析命令行参数
    mode = "all"  # 默认全RNA模式
    similarity_threshold = 70  # 兜底默认值（工作流已改为40%）
    
    if len(sys.argv) > 1:
        mode = sys.argv[1].lower()
    if len(sys.argv) > 2:
        similarity_threshold = int(sys.argv[2])
    
    print("="*70)
    if mode == "antibacterial":
        print("🦠 全自动新型抗菌药物靶点发现与设计引擎")
        targets_df = load_antibacterial_targets()
        output_dir = ANTIBACTERIAL_OUTPUT_DIR
        report_title = "新型抗菌药物靶点发现与设计报告"
        report_subtitle = "基于核糖开关和核糖体RNA结构的AIDD药物重定位分析"
    else:
        print("🧬 全RNA靶点药物重定位自动化分析引擎")
        targets_df = load_all_rna_targets()
        output_dir = ALL_RNA_OUTPUT_DIR
        report_title = "全RNA靶点药物重定位分析报告"
        report_subtitle = "基于所有PDB RNA结构的系统性药物重定位分析"
    print("="*70)
    
    if len(targets_df) == 0:
        print("❌ 未找到任何有效靶点")
        return
    
    # 批量分析
    result_df = batch_analyze_targets(
        targets_df, 
        similarity_threshold, 
        output_dir, 
        report_title, 
        report_subtitle
    )
    
    print("\n" + "="*70)
    print("✅ 分析流程完成！")
    print("="*70)
    print(f"\n📂 所有结果已保存到 '{output_dir}' 目录")

if __name__ == "__main__":
    main()
