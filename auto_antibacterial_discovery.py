import pandas as pd
import requests
import urllib.parse
import time
import json
import os
from functools import lru_cache
from datetime import datetime

# =============================================================================
# 🦠 全自动新型抗菌药物靶点发现与设计引擎
# 基于RNA结构与AIDD药物重定位技术
# 自动完成：靶点筛选 → 批量药物分析 → 优先级排序 → 生成专业报告
# =============================================================================

# --------------------------
# 配置与缓存系统（与原平台完全一致）
# --------------------------
SMILES_CACHE_FILE = "smiles_cache.json"
DRUG_SEARCH_CACHE_FILE = "drug_search_cache.json"
OUTPUT_DIR = "antibacterial_results"

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ATC分类与适应症映射（与原平台完全一致）
ATC_FULL_MAP = {
    "J01": "全身用抗菌药", "J02": "全身用抗真菌药", "J04": "抗分枝杆菌药", "J05": "全身用抗病毒药",
    "A07": "肠道抗感染药和肠道消炎药"
}

INDICATION_CN_MAP = {
    "bacterial infection": "细菌感染", "infection": "感染性疾病", "tuberculosis": "结核病",
    "pneumonia": "肺炎", "sepsis": "脓毒症", "urinary tract infection": "尿路感染"
}

MANUAL_DRUG_INDICATION = {
    "PAROMOMYCIN": "肠道阿米巴病、细菌性痢疾",
    "PAROMOMYCIN SULFATE": "肠道阿米巴病、细菌性痢疾",
    "NETILMICIN": "敏感菌所致的呼吸道、泌尿道、皮肤软组织感染",
    "NETILMICIN SULFATE": "敏感菌所致的呼吸道、泌尿道、皮肤软组织感染",
    "KANAMYCIN": "敏感菌所致的严重感染",
    "KANAMYCIN SULFATE": "敏感菌所致的严重感染",
    "GENTAMICIN": "革兰氏阴性菌所致的严重感染",
    "GENTAMICIN SULFATE": "革兰氏阴性菌所致的严重感染",
    "TOBRAMYCIN": "铜绿假单胞菌等革兰氏阴性菌感染",
    "TOBRAMYCIN SULFATE": "铜绿假单胞菌等革兰氏阴性菌感染",
    "AMIKACIN": "敏感菌所致的严重感染",
    "AMIKACIN SULFATE": "敏感菌所致的严重感染",
    "STREPTOMYCIN": "结核病、鼠疫",
    "STREPTOMYCIN SULFATE": "结核病、鼠疫",
    "ERYTHROMYCIN": "呼吸道感染、皮肤软组织感染",
    "TETRACYCLINE": "立克次体病、支原体肺炎"
}

# --------------------------
# 通用工具函数
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
# 核心函数（与原平台完全一致，保证结果一致性）
# --------------------------
@lru_cache(maxsize=1000)
def get_smiles_by_id(ligand_id):
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
                
                # 智能适应症获取
                indication_list = []
                if drug_name in MANUAL_DRUG_INDICATION:
                    indication_list.append(MANUAL_DRUG_INDICATION[drug_name])
                
                if not indication_list:
                    drug_indications = m.get('drug_indications', [])
                    for ind in drug_indications:
                        mesh_term = ind.get('mesh_heading', '').lower()
                        efo_term = ind.get('efo_term', '').lower()
                        for eng, cn in INDICATION_CN_MAP.items():
                            if eng in mesh_term or eng in efo_term:
                                indication_list.append(cn)
                                break
                
                if not indication_list:
                    base_name = drug_name.replace(' SULFATE', '').replace(' HYDROCHLORIDE', '')
                    if base_name in MANUAL_DRUG_INDICATION:
                        indication_list.append(MANUAL_DRUG_INDICATION[base_name])
                
                if not indication_list:
                    indication_list.append("细菌感染")
                
                # 药物分类
                drug_class = "全身用抗菌药"
                atc_list = m.get('atc_classifications', [])
                if atc_list:
                    atc_code = atc_list[0][:3]
                    if atc_code in ATC_FULL_MAP:
                        drug_class = ATC_FULL_MAP[atc_code]
                
                drugs.append({
                    "药物名称": m.get('pref_name'),
                    "ChEMBL ID": m.get('molecule_chembl_id'),
                    "相似度(%)": round(float(m.get('similarity', 0)), 2),
                    "分子量": round(float(m.get('molecule_properties', {}).get('full_mwt', 0)), 2),
                    "药物分类": drug_class,
                    "治疗适应症": "、".join(list(set(indication_list)))
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
# 自动化抗菌靶点筛选引擎
# --------------------------
def load_and_filter_antibacterial_targets():
    """自动加载并筛选抗菌RNA靶点（核糖开关和核糖体）"""
    print("🔍 正在加载并筛选抗菌RNA靶点...")
    
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
# 批量药物分析与优先级排序
# --------------------------
def batch_analyze_antibacterial_targets(targets_df, similarity_threshold=70):
    """批量分析所有抗菌靶点的潜在药物"""
    print("\n🚀 开始批量药物重定位分析...")
    
    all_results = []
    total_targets = len(targets_df)
    
    for idx, (_, row) in enumerate(targets_df.iterrows()):
        pdb_id = row['PDB ID']
        ligand_id = row['MainLigandID']
        category = row['Category']
        description = row['Description (描述)']
        
        print(f"   正在处理 [{idx+1}/{total_targets}] {pdb_id} (配体: {ligand_id})")
        
        smiles = get_smiles_by_id(ligand_id)
        if not smiles:
            print(f"   ⚠️ {pdb_id} 无法获取SMILES，已跳过")
            continue
        
        drugs, msg, total = search_chembl_drugs(smiles, similarity_threshold)
        
        for drug in drugs:
            drug["RNA类别"] = category
            drug["PDB ID"] = pdb_id
            drug["配体ID"] = ligand_id
            drug["RNA结构描述"] = description
            all_results.append(drug)
    
    print(f"\n✅ 批量分析完成！共匹配到 {len(all_results)} 条抗菌药物记录")
    
    # 优先级排序
    print("\n📊 正在进行药物优先级排序...")
    result_df = pd.DataFrame(all_results)
    
    # 计算优先级分数
    def calculate_priority(row):
        score = 0
        # 相似度权重：70%
        score += row["相似度(%)"] * 0.7
        # 药物分类权重：20%
        if row["药物分类"] == "全身用抗菌药":
            score += 20
        # 细菌相关靶点权重：10%
        if "bacterial" in str(row["RNA结构描述"]).lower():
            score += 10
        return round(score, 2)
    
    result_df["优先级分数"] = result_df.apply(calculate_priority, axis=1)
    result_df = result_df.sort_values(by=["优先级分数", "相似度(%)"], ascending=False)
    
    # 重新排列列顺序
    column_order = [
        "优先级分数", "RNA类别", "PDB ID", "配体ID", "药物名称", 
        "药物分类", "治疗适应症", "相似度(%)", "分子量", "ChEMBL ID", "RNA结构描述"
    ]
    result_df = result_df[column_order]
    
    return result_df

# --------------------------
# 生成专业HTML报告
# --------------------------
def generate_html_report(result_df, targets_df):
    """生成完整的抗菌药物发现报告"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # 统计信息
    total_targets = len(targets_df)
    total_drugs = len(result_df)
    unique_drugs = result_df["药物名称"].nunique()
    top_10_drugs = result_df.head(10).to_html(index=False, classes="table table-striped table-hover")
    
    # 按RNA类别统计
    riboswitch_drugs = len(result_df[result_df["RNA类别"] == "核糖开关 (Riboswitch)"])
    rrna_drugs = len(result_df[result_df["RNA类别"] == "核糖体 (rRNA)"])
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="zh-CN">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>新型抗菌药物靶点发现与设计报告</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; }}
            .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2rem 0; margin-bottom: 2rem; }}
            .section {{ margin-bottom: 2rem; padding: 1.5rem; background: #f8f9fa; border-radius: 10px; }}
            .stat-card {{ background: white; padding: 1.5rem; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); text-align: center; }}
            .stat-number {{ font-size: 2rem; font-weight: bold; color: #667eea; }}
            table {{ font-size: 0.9rem; }}
            .footer {{ text-align: center; padding: 2rem; color: #6c757d; border-top: 1px solid #eee; margin-top: 3rem; }}
        </style>
    </head>
    <body>
        <div class="header">
            <div class="container">
                <h1 class="text-center">🦠 新型抗菌药物靶点发现与设计报告</h1>
                <p class="text-center lead">基于RNA结构与AIDD药物重定位技术</p>
                <p class="text-center">报告生成时间: {timestamp}</p>
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
                            <div class="stat-number">{riboswitch_drugs}/{rrna_drugs}</div>
                            <div>核糖开关/核糖体药物</div>
                        </div>
                    </div>
                </div>
            </div>

            <div class="section">
                <h2>🏆 高优先级候选药物 TOP 10</h2>
                <p>基于相似度、药物分类和靶点相关性综合排序</p>
                {top_10_drugs}
            </div>

            <div class="section">
                <h2>💡 实验验证建议</h2>
                <ol>
                    <li><strong>最低抑菌浓度(MIC)测定</strong>：首先测试TOP 10药物对目标细菌的体外抗菌活性</li>
                    <li><strong>结合亲和力验证</strong>：使用等温滴定量热法(ITC)或表面等离子体共振(SPR)验证药物与RNA的结合</li>
                    <li><strong>作用机制研究</strong>：通过报告基因实验验证药物对RNA功能的调控作用</li>
                    <li><strong>动物体内实验</strong>：在小鼠感染模型中验证药物的体内疗效</li>
                    <li><strong>耐药性评估</strong>：测试药物诱导细菌产生耐药性的难易程度</li>
                </ol>
            </div>

            <div class="section">
                <h2>🚀 后续研究方向</h2>
                <ul>
                    <li><strong>老药新用</strong>：优先开发已上市药物的新适应症，缩短研发周期</li>
                    <li><strong>结构优化</strong>：基于RNA-配体复合物结构进行药物分子修饰，提高亲和力和特异性</li>
                    <li><strong>联合用药</strong>：与现有抗生素联合使用，逆转细菌耐药性</li>
                    <li><strong>兽用抗生素开发</strong>：针对养殖业常见致病菌开发专用药物</li>
                </ul>
            </div>
        </div>

        <div class="footer">
            <p>本报告由RNA结构精细化分类与AIDD药物重定位平台自动生成</p>
            <p>© 2025 抗菌药物研发自动化系统</p>
        </div>
    </body>
    </html>
    """
    
    report_path = os.path.join(OUTPUT_DIR, "antibacterial_drug_discovery_report.html")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    csv_path = os.path.join(OUTPUT_DIR, "antibacterial_candidate_drugs.csv")
    result_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    
    print(f"\n📄 报告生成完成！")
    print(f"   - 完整报告: {report_path}")
    print(f"   - 候选药物列表: {csv_path}")
    
    return report_path, csv_path

# --------------------------
# 主函数：一键运行全流程
# --------------------------
def main():
    print("="*60)
    print("🦠 全自动新型抗菌药物靶点发现与设计引擎")
    print("="*60)
    
    # 步骤1：筛选抗菌靶点
    targets_df = load_and_filter_antibacterial_targets()
    
    if len(targets_df) == 0:
        print("❌ 未找到任何抗菌RNA靶点")
        return
    
    # 步骤2：批量药物分析
    result_df = batch_analyze_antibacterial_targets(targets_df, similarity_threshold=70)
    
    if len(result_df) == 0:
        print("❌ 未匹配到任何抗菌药物")
        return
    
    # 步骤3：生成报告
    report_path, csv_path = generate_html_report(result_df, targets_df)
    
    print("\n" + "="*60)
    print("✅ 全自动抗菌药物发现流程完成！")
    print("="*60)
    print(f"\n📌 高优先级候选药物TOP 3:")
    for i in range(min(3, len(result_df))):
        row = result_df.iloc[i]
        print(f"   {i+1}. {row['药物名称']} (相似度: {row['相似度(%)']}%, 靶点: {row['PDB ID']})")
    
    print(f"\n📂 所有结果已保存到 '{OUTPUT_DIR}' 目录")
    print(f"🌐 请打开 {report_path} 查看完整报告")

if __name__ == "__main__":
    main()
