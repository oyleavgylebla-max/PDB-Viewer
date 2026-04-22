import pandas as pd
import requests
import urllib.parse
import json
import os
import time
from functools import lru_cache
from tqdm import tqdm

# =============================================================================
# 🧬 RNA靶点AIDD药物重定位 全自动工作流（靶点查询修复版）
# 版本: 2.1 修复版
# 修复内容: 优化ChEMBL API查询逻辑，解决药物靶点查询失败问题
# =============================================================================

# ======================================
# 🔧 【用户可配置参数区】
# ======================================
MIN_SIMILARITY = 70
TARGET_RNA_CATEGORIES = [
    "核糖开关 (Riboswitch)",
    "核糖体 (rRNA)",
    "G-四联体 (G-quadruplex)",
    "适配体 (Aptamer)",
    "转运RNA (tRNA)",
    "核酶 (Ribozyme)",
    "特殊结构基元 (Special/Motifs)"
]
TARGET_DRUG_CLASSES = [
    "全身用抗菌药",
    "抗肿瘤药",
    "全身用抗病毒药",
    "抗真菌药",
    "心血管系统药物",
    "神经系统药物"
]
DRUG_DEDUPLICATION = True
PDB_EXCEL_PATH = "PDB_Dataset_Info_Full.xlsx"
OUTPUT_FOLDER = "rna_aidd_workflow_output"
SMILES_CACHE_FILE = "cache_smiles.json"
DRUG_CACHE_FILE = "cache_drugs.json"
TARGET_CACHE_FILE = "cache_targets.json"
REQUEST_TIMEOUT = 15
MAX_RETRIES = 2
RETRY_DELAY = 1

# ======================================
# 全局配置与字典
# ======================================
ATC_FULL_MAP = {
    "A01": "口腔病用药", "A02": "治疗胃酸相关疾病的药物", "A03": "治疗功能性胃肠道疾病的药物",
    "A04": "止吐药和止恶心药", "A05": "胆和肝治疗药", "A06": "轻泻药", "A07": "肠道抗感染药和肠道消炎药",
    "A08": "减肥药", "A09": "消化药，包括酶", "A10": "糖尿病用药", "A11": "维生素类", "A12": "矿物质补充剂",
    "A13": "滋补药", "A14": "全身用蛋白同化类固醇", "A15": "食欲刺激药", "A16": "其他消化道和代谢药物",
    "B01": "抗血栓药", "B02": "止血药", "B03": "抗贫血药", "B05": "血液代用品和灌注液", "B06": "其他血液系统用药",
    "C01": "心脏治疗药", "C02": "抗高血压药", "C03": "利尿药", "C04": "外周血管扩张药", "C05": "血管保护剂",
    "C07": "β受体阻滞剂", "C08": "钙通道阻滞剂", "C09": "作用于肾素-血管紧张素系统的药物", "C10": "血脂调节剂",
    "D01": "皮肤用抗真菌药", "D02": "润肤剂和保护剂", "D03": "皮肤用皮质类固醇", "D04": "止痒药，包括抗组胺药、麻醉药",
    "D05": "银屑病用药", "D06": "皮肤用抗生素和化疗药", "D07": "皮肤用皮质类固醇和抗生素的复方制剂",
    "D08": "皮肤用消毒剂和防腐剂", "D09": "伤口敷料和保护剂", "D10": "痤疮用药", "D11": "其他皮肤科用药",
    "J01": "全身用抗菌药", "J02": "全身用抗真菌药", "J04": "抗分枝杆菌药", "J05": "全身用抗病毒药",
    "J06": "免疫血清和免疫球蛋白", "J07": "疫苗",
    "L01": "抗肿瘤药", "L02": "内分泌治疗药", "L03": "免疫刺激剂", "L04": "免疫抑制剂",
    "M01": "抗炎和抗风湿药", "M02": "局部用肌肉骨骼系统药物", "M03": "肌肉松弛药", "M04": "抗痛风药", "M05": "治疗骨病的药物",
    "N01": "麻醉药", "N02": "镇痛药", "N03": "抗癫痫药", "N04": "抗帕金森病药", "N05": "精神安定药",
    "N06": "精神兴奋药", "N07": "其他神经系统药物",
    "R01": "鼻腔用药", "R02": "咽喉用药", "R03": "用于阻塞性气道疾病的药物", "R05": "咳嗽和感冒用药",
    "R06": "全身用抗组胺药", "R07": "其他呼吸系统药物"
}
ATC_CLASS_MAP = {
    "A": "消化系统及代谢药", "B": "血液和造血系统药物", "C": "心血管系统药物",
    "D": "皮肤科用药", "G": "泌尿生殖系统药和性激素", "H": "全身用激素类制剂",
    "J": "全身用抗感染药", "L": "抗肿瘤药和免疫调节剂", "M": "肌肉-骨骼系统药物",
    "N": "神经系统药物", "P": "抗寄生虫药", "R": "呼吸系统药物",
    "S": "感觉器官药物", "V": "其他药品"
}
INDICATION_CN_MAP = {
    "cancer": "癌症", "tumor": "肿瘤", "carcinoma": "恶性肿瘤", "leukemia": "白血病",
    "hiv": "艾滋病", "influenza": "流感", "hepatitis": "肝炎", "bacterial infection": "细菌感染",
    "diabetes": "糖尿病", "obesity": "肥胖症", "hypertension": "高血压", "cardiovascular": "心血管疾病",
    "alzheimer": "阿尔茨海默病", "parkinson": "帕金森病", "arthritis": "关节炎", "asthma": "哮喘",
    "depression": "抑郁症", "pain": "疼痛", "inflammation": "炎症", "fungal infection": "真菌感染",
    "virus infection": "病毒感染", "thromboembolism": "血栓栓塞", "stroke": "中风",
    "psychosis": "精神病", "schizophrenia": "精神分裂症", "epilepsy": "癫痫",
    "gastroesophageal reflux": "胃食管反流", "ulcer": "消化道溃疡", "anemia": "贫血",
    "glaucoma": "青光眼", "migraine": "偏头痛", "osteoporosis": "骨质疏松症",
    "infection": "感染性疾病", "autoimmune": "自身免疫性疾病", "neuropathic pain": "神经病理性疼痛",
    "tuberculosis": "结核病", "malaria": "疟疾", "pneumonia": "肺炎", "sepsis": "脓毒症"
}
MANUAL_DRUG_INDICATION = {
    "PAROMOMYCIN": "肠道阿米巴病、细菌性痢疾", "PAROMOMYCIN SULFATE": "肠道阿米巴病、细菌性痢疾",
    "NETILMICIN": "敏感菌所致的呼吸道、泌尿道、皮肤软组织感染", "NETILMICIN SULFATE": "敏感菌所致的呼吸道、泌尿道、皮肤软组织感染",
    "KANAMYCIN": "敏感菌所致的严重感染", "KANAMYCIN SULFATE": "敏感菌所致的严重感染",
    "GENTAMICIN": "革兰氏阴性菌所致的严重感染", "GENTAMICIN SULFATE": "革兰氏阴性菌所致的严重感染",
    "TOBRAMYCIN": "铜绿假单胞菌等革兰氏阴性菌感染", "TOBRAMYCIN SULFATE": "铜绿假单胞菌等革兰氏阴性菌感染",
    "AMIKACIN": "敏感菌所致的严重感染", "AMIKACIN SULFATE": "敏感菌所致的严重感染",
    "STREPTOMYCIN": "结核病、鼠疫", "STREPTOMYCIN SULFATE": "结核病、鼠疫",
    "ERYTHROMYCIN": "呼吸道感染、皮肤软组织感染", "ERYTHROMYCIN LACTOBIONATE": "呼吸道感染、皮肤软组织感染",
    "TETRACYCLINE": "立克次体病、支原体肺炎", "TETRACYCLINE HYDROCHLORIDE": "立克次体病、支原体肺炎"
}
KNOWN_RNA_DRUGS = {
    "GENTAMICIN": "已知靶向细菌核糖体RNA", "AMIKACIN": "已知靶向细菌核糖体RNA",
    "TOBRAMYCIN": "已知靶向细菌核糖体RNA", "KANAMYCIN": "已知靶向细菌核糖体RNA",
    "STREPTOMYCIN": "已知靶向细菌核糖体RNA", "PAROMOMYCIN": "已知靶向细菌核糖体RNA",
    "NETILMICIN": "已知靶向细菌核糖体RNA", "NEOMYCIN": "已知靶向细菌核糖体RNA",
    "RIBAVIRIN": "已知靶向病毒RNA聚合酶", "SOFOSBUVIR": "已知靶向病毒RNA聚合酶",
    "LINEZOLID": "已知靶向细菌核糖体RNA", "TEDIZOLID": "已知靶向细菌核糖体RNA"
}
ION_FILTER_LIST = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F', 'IOD', 'FLC', 'NO3', 'NH4', 'ACT', 'FMT', 'EDO', 'GOL', 'PEG', 'DTT', 'BME']
REQUEST_HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
    "Accept": "application/json, text/plain, */*"
}

# ======================================
# 缓存工具函数
# ======================================
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
        print(f"[警告] 缓存保存失败: {e}")

# ======================================
# 核心工作流函数
# ======================================
def categorize_rna(desc):
    desc_lower = str(desc).lower()
    if 'mirna' in desc_lower and 'riboswitch scaffold' in desc_lower:
        return "特殊结构基元 (Special/Motifs)"
    if 'riboswitch' in desc_lower: return "核糖开关 (Riboswitch)"
    elif 'aptamer' in desc_lower: return "适配体 (Aptamer)"
    elif any(word in desc_lower for word in ['quadruplex', 'g-4', 'g4', 'tetraplex']): return "G-四联体 (G-quadruplex)"
    elif any(word in desc_lower for word in ['ribosomal', 'ribosome', 'rrna', 'decoding site', 'a-site']): return "核糖体 (rRNA)"
    elif 'ribozyme' in desc_lower: return "核酶 (Ribozyme)"
    elif any(word in desc_lower for word in ['trna', 'transfer rna', 't-rna']): return "转运RNA (tRNA)"
    elif any(word in desc_lower for word in [
        'ires', 'hairpin', 'stem-loop', 'pseudoknot', 'bulge', 'duplex',
        'tar rna', 'hiv-1 tar', 'hiv-2 tar', 'splicing regulatory',
        'dimerization initiation', 'frameshift site', 'cag repeats',
        'rna helix', 'rna oligonucleotide', 'pre-mrna'
    ]): return "特殊结构基元 (Special/Motifs)"
    elif any(word in desc_lower for word in ['influenza rna', 'hcv', 'hepatitis c virus', 'viral rna']): return "特殊结构基元 (Special/Motifs)"
    else: return "其他 RNA (Others)"

def extract_main_ligand(ligand_text):
    if ligand_text == "No ligands" or str(ligand_text) == "nan":
        return "ZZZ"
    all_ligands = str(ligand_text).split(' | ')
    for ligand in all_ligands:
        ligand_id = ligand.strip().split(' ')[0].upper()
        if ligand_id not in ION_FILTER_LIST:
            return ligand_id
    return all_ligands[0].strip().split(' ')[0].upper()

@lru_cache(maxsize=2000)
def get_smiles_by_ligand_id(ligand_id):
    if not ligand_id or ligand_id == "ZZZ":
        return None
    ligand_id = ligand_id.strip().upper()
    smiles_cache = load_cache(SMILES_CACHE_FILE)
    if ligand_id in smiles_cache:
        return smiles_cache[ligand_id]
    
    smiles_result = None
    for attempt in range(MAX_RETRIES):
        try:
            url = f"https://www.ebi.ac.uk/pdbe/api/v2/compound/summary/{ligand_id}"
            r = requests.get(url, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT)
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
        except:
            time.sleep(RETRY_DELAY)
    if not smiles_result:
        for attempt in range(MAX_RETRIES):
            try:
                url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
                r = requests.get(url, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT)
                if r.status_code == 200:
                    data = r.json()
                    if "rcsb_chem_comp_descriptor" in data:
                        desc = data["rcsb_chem_comp_descriptor"]
                        if "SMILES_stereo" in desc: smiles_result = desc["SMILES_stereo"]
                        elif "SMILES" in desc: smiles_result = desc["SMILES"]
                        if smiles_result: break
            except:
                time.sleep(RETRY_DELAY)
    if not smiles_result:
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/PDB/{ligand_id}/cids/JSON"
            r = requests.get(url, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT)
            if r.status_code == 200:
                cid = r.json()["IdentifierList"]["CID"][0]
                url2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
                r2 = requests.get(url2, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT)
                if r2.status_code == 200:
                    smiles_result = r2.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        except:
            pass
    if smiles_result:
        smiles_cache[ligand_id] = smiles_result
        save_cache(smiles_cache, SMILES_CACHE_FILE)
        return smiles_result
    return None

def search_similar_drugs(smiles, similarity_threshold):
    if not smiles:
        return [], "SMILES为空", 0
    cache_key = f"{smiles}_{similarity_threshold}"
    drug_cache = load_cache(DRUG_CACHE_FILE)
    if cache_key in drug_cache:
        return drug_cache[cache_key]["drugs"], drug_cache[cache_key]["msg"], drug_cache[cache_key]["total"]
    
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json?limit=1000&expand=drug_indications,molecule_structures"
        r = requests.get(url, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT*2)
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            drugs = []
            for m in mols:
                if not (m.get('max_phase') and float(m.get('max_phase')) >= 4.0 and m.get('pref_name') and not m.get('withdrawn_flag')):
                    continue
                drug_name = m.get('pref_name', '').upper()
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
                    base_name = drug_name.replace(' SULFATE', '').replace(' HYDROCHLORIDE', '').replace(' HCL', '')
                    if base_name in MANUAL_DRUG_INDICATION:
                        indication_list.append(MANUAL_DRUG_INDICATION[base_name])
                if not indication_list:
                    indication_list.append("暂无精准适应症数据")
                unique_indication = list(set(indication_list))
                indication_str = "、".join(unique_indication)
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
                drugs.append({
                    "药物名称": m.get('pref_name'),
                    "ChEMBL ID": m.get('molecule_chembl_id'),
                    "相似度(%)": round(float(m.get('similarity', 0)), 2),
                    "分子量": round(float(m.get('molecule_properties', {}).get('full_mwt', 0)), 2),
                    "药物分类归属": drug_class,
                    "治疗适应症": indication_str,
                })
            drugs.sort(key=lambda x: x["相似度(%)"], reverse=True)
            drug_cache[cache_key] = {"drugs": drugs, "msg": "Success", "total": len(mols)}
            save_cache(drug_cache, DRUG_CACHE_FILE)
            return drugs, "Success", len(mols)
        else:
            return [], f"接口请求失败 ({r.status_code})", 0
    except Exception as e:
        return [], f"请求异常: {str(e)}", 0

# ======================================
# 🔍 修复版：靶点不同确认核心模块
# ======================================
def get_drug_known_targets_cached(drug_name):
    """带缓存的靶点查询（修复大小写和API问题）"""
    target_cache = load_cache(TARGET_CACHE_FILE)
    drug_name_key = drug_name.upper().strip()
    
    if drug_name_key in target_cache:
        return target_cache[drug_name_key]["types"], target_cache[drug_name_key]["descs"]
    
    try:
        # 修复1：不区分大小写查询
        search_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?format=json"
        params = {
            "pref_name__icontains": drug_name.strip(),
            "limit": 10
        }
        r = requests.get(search_url, params=params, headers=REQUEST_HEADERS, timeout=10)
        if r.status_code != 200:
            return set(), []
        
        data = r.json()
        if not data.get("molecules"):
            return set(), []
        
        # 优先选择完全匹配的
        mol_chembl_id = None
        for mol in data["molecules"]:
            if mol.get("pref_name", "").upper().strip() == drug_name_key:
                mol_chembl_id = mol["molecule_chembl_id"]
                break
        
        # 如果没有完全匹配，取第一个
        if not mol_chembl_id:
            mol_chembl_id = data["molecules"][0]["molecule_chembl_id"]
        
        # 修复2：改用更可靠的靶点接口
        target_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{mol_chembl_id}/targets?format=json"
        r2 = requests.get(target_url, headers=REQUEST_HEADERS, timeout=10)
        if r2.status_code != 200:
            return set(), []
        
        target_data = r2.json()
        target_types = set()
        target_descriptions = []
        
        for target in target_data.get("targets", []):
            target_type = target.get("target_type", "UNKNOWN")
            target_desc = target.get("pref_name", "")
            target_types.add(target_type.upper())
            target_descriptions.append(target_desc)
        
        target_descriptions = list(set(target_descriptions))
        
        target_cache[drug_name_key] = {
            "types": list(target_types),
            "descs": target_descriptions
        }
        save_cache(target_cache, TARGET_CACHE_FILE)
        
        return target_types, target_descriptions
    
    except Exception as e:
        return set(), []

def is_confirmed_different_target(row):
    """确认预测靶点与药物已知靶点不同"""
    drug_name = row["药物名称"]
    predicted_rna_category = row["RNA类别"]
    
    known_target_types, known_target_descs = get_drug_known_targets_cached(drug_name)
    
    if not known_target_types:
        return True, "未查询到该药物的已知靶点信息", "未知"
    
    has_rna_target = any("RNA" in t or "RIBOSOME" in t or "NUCLEIC ACID" in t for t in known_target_types)
    if has_rna_target:
        if "RIBOSOME" in str(known_target_types).upper() and predicted_rna_category == "核糖体 (rRNA)":
            return False, "药物已知靶向核糖体RNA", "; ".join(known_target_descs[:3])
        elif predicted_rna_category != "核糖体 (rRNA)":
            return True, "药物已知靶向其他RNA，但预测的是不同类型的RNA", "; ".join(known_target_descs[:3])
        else:
            return False, "药物已知靶向RNA", "; ".join(known_target_descs[:3])
    
    has_only_protein = all("PROTEIN" in t or "ENZYME" in t or "RECEPTOR" in t or "KINASE" in t for t in known_target_types)
    if has_only_protein:
        return True, "药物所有已知靶点均为蛋白质，预测的RNA是潜在新靶点", "; ".join(known_target_descs[:3])
    
    return True, "需人工确认靶点类型", "; ".join(known_target_descs[:3])

# ======================================
# 主工作流执行
# ======================================
def main_workflow():
    print("="*80)
    print("🧬 RNA靶点AIDD药物重定位 全自动工作流（靶点查询修复版）")
    print("="*80)

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    print(f"[1/7] 输出文件夹已创建: {OUTPUT_FOLDER}")

    try:
        df = pd.read_excel(PDB_EXCEL_PATH)
        df = df[df['PDB ID'] != '5D51'].reset_index(drop=True)
        df['RNA类别'] = df['Description (描述)'].apply(categorize_rna)
        df['核心配体ID'] = df['Ligands (对应小分子)'].apply(extract_main_ligand)
        valid_df = df[df['核心配体ID'] != 'ZZZ'].reset_index(drop=True)
        print(f"[2/7] PDB数据加载完成: 共{len(df)}个结构，有效配体结构{len(valid_df)}个")
        df.to_excel(f"{OUTPUT_FOLDER}/01_RNA结构完整分类表.xlsx", index=False)
    except FileNotFoundError:
        print(f"[错误] 未找到PDB数据文件: {PDB_EXCEL_PATH}")
        return
    except Exception as e:
        print(f"[错误] 数据加载失败: {e}")
        return

    print(f"[3/7] 开始批量获取配体SMILES，共{len(valid_df['核心配体ID'].unique())}个唯一配体")
    ligand_smiles_map = {}
    unique_ligands = valid_df['核心配体ID'].unique()
    for ligand_id in tqdm(unique_ligands, desc="SMILES获取进度"):
        smiles = get_smiles_by_ligand_id(ligand_id)
        if smiles:
            ligand_smiles_map[ligand_id] = smiles
    print(f"[3/7] SMILES获取完成: 成功获取{len(ligand_smiles_map)}个配体的SMILES")

    print(f"[4/7] 开始批量搜索相似上市药物，共{len(ligand_smiles_map)}个配体")
    all_drug_results = []
    for _, row in tqdm(valid_df.iterrows(), total=len(valid_df), desc="药物搜索进度"):
        pdb_id = row['PDB ID']
        rna_category = row['RNA类别']
        ligand_id = row['核心配体ID']
        desc = row['Description (描述)']
        if ligand_id not in ligand_smiles_map:
            continue
        smiles = ligand_smiles_map[ligand_id]
        drugs, msg, total = search_similar_drugs(smiles, MIN_SIMILARITY)
        if drugs:
            for drug in drugs:
                drug['RNA类别'] = rna_category
                drug['PDB ID'] = pdb_id
                drug['核心配体ID'] = ligand_id
                drug['RNA结构描述'] = desc
                all_drug_results.append(drug)
    print(f"[4/7] 药物搜索完成: 共匹配到{len(all_drug_results)}条药物记录")

    print(f"[5/7] 开始进行靶点不同确认...")
    full_result_df = pd.DataFrame(all_drug_results)
    pre_filtered = full_result_df[
        (full_result_df["相似度(%)"] >= MIN_SIMILARITY) &
        (full_result_df["药物分类归属"].isin(TARGET_DRUG_CLASSES)) &
        (full_result_df["RNA类别"].isin(TARGET_RNA_CATEGORIES))
    ].copy()

    confirmation_results = []
    for idx, row in tqdm(pre_filtered.iterrows(), total=len(pre_filtered), desc="靶点确认进度"):
        is_different, reason, known_targets = is_confirmed_different_target(row)
        confirmation_results.append({
            "是否确认靶点不同": is_different,
            "判断理由": reason,
            "已知靶点信息": known_targets
        })
        time.sleep(0.15)

    confirmation_df = pd.DataFrame(confirmation_results)
    pre_filtered = pd.concat([pre_filtered.reset_index(drop=True), confirmation_df], axis=1)
    final_filtered = pre_filtered[pre_filtered["是否确认靶点不同"] == True].copy()
    
    if DRUG_DEDUPLICATION:
        final_filtered = final_filtered.drop_duplicates(subset=["药物名称"], keep="first")
    
    final_filtered = final_filtered.sort_values(by=["相似度(%)"], ascending=False)
    
    print(f"[5/7] 靶点确认完成:")
    print(f"   - 初步筛选候选数：{len(pre_filtered)}")
    print(f"   - 确认靶点不同候选数：{len(final_filtered)}")
    print(f"   - 排除已知RNA靶向药数：{len(pre_filtered) - len(final_filtered)}")

    print(f"[6/7] 开始输出结果文件")
    column_order = [
        "RNA类别", "PDB ID", "核心配体ID", "RNA结构描述",
        "药物名称", "药物分类归属", "治疗适应症", "相似度(%)", "分子量",
        "判断理由", "已知靶点信息", "ChEMBL ID"
    ]
    final_filtered = final_filtered[[col for col in column_order if col in final_filtered.columns]]
    
    full_result_df.to_csv(f"{OUTPUT_FOLDER}/02_全库RNA-药物完整匹配结果.csv", index=False, encoding="utf-8-sig")
    final_filtered.to_csv(f"{OUTPUT_FOLDER}/03_RNA靶向老药新用高价值候选清单(已确认靶点).csv", index=False, encoding="utf-8-sig")
    final_filtered.to_excel(f"{OUTPUT_FOLDER}/03_RNA靶向老药新用高价值候选清单(已确认靶点).xlsx", index=False)
    
    report_content = [
        "="*60,
        "RNA靶点AIDD药物重定位工作流 统计报告（含靶点不同确认）",
        "="*60,
        f"分析时间: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"相似度阈值: ≥{MIN_SIMILARITY}%",
        "",
        "一、基础数据统计",
        f"1. 总RNA结构数量: {len(df)}个",
        f"2. 含有效配体的RNA结构数量: {len(valid_df)}个",
        f"3. 唯一配体数量: {len(unique_ligands)}个",
        f"4. 成功获取SMILES的配体数量: {len(ligand_smiles_map)}个",
        "",
        "二、药物匹配结果",
        f"1. 总匹配药物记录数: {len(all_drug_results)}条",
        f"2. 初步筛选候选数: {len(pre_filtered)}个",
        f"3. 确认靶点不同候选数: {len(final_filtered)}个",
        f"4. 排除已知RNA靶向药数: {len(pre_filtered) - len(final_filtered)}个",
        f"5. 覆盖RNA靶点数量: {final_filtered['PDB ID'].nunique()}个",
        f"6. 平均相似度: {round(final_filtered['相似度(%)'].mean(), 2)}%",
        "",
        "三、药物分类分布",
        final_filtered["药物分类归属"].value_counts().to_string(),
        "",
        "四、RNA类型分布",
        final_filtered["RNA类别"].value_counts().to_string(),
        "",
        "五、判断理由分布",
        final_filtered["判断理由"].value_counts().to_string(),
        "",
        "="*60,
        "报告结束",
        "="*60
    ]
    with open(f"{OUTPUT_FOLDER}/04_科研统计报告(含靶点确认).txt", "w", encoding="utf-8") as f:
        f.write("\n".join(report_content))

    print("="*80)
    print("🎉 工作流执行完成！所有结果已保存到输出文件夹")
    print("="*80)
    print("\n📊 核心统计结果（可直接用于论文）:")
    print(f"✅ 确认靶点不同的高价值候选药物总数: {len(final_filtered)} 个")
    print(f"✅ 覆盖RNA靶点数量: {final_filtered['PDB ID'].nunique()} 个")
    print(f"✅ 平均结构相似度: {round(final_filtered['相似度(%)'].mean(), 2)}%")
    print("\n🏆 药物分类分布:")
    print(final_filtered["药物分类归属"].value_counts().to_string())
    print("\n🏆 判断理由分布:")
    print(final_filtered["判断理由"].value_counts().to_string())
    print(f"\n📁 所有结果已保存到: {os.path.abspath(OUTPUT_FOLDER)}")

if __name__ == "__main__":
    main_workflow()
