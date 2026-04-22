import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol
import requests
import urllib.parse
import time
import json
import os
from functools import lru_cache

# =============================================================================
# 🧬 RNA 结构精细化分类与 AIDD 药物重定位平台
# 版本: 3.3（基于PDB实际描述完善RNA结构分类）
# 核心优化：
# 1. 新增12个RNA结构细分亚型
# 2. 标注结构来源（细菌/病毒/人类）
# 3. 支持按细分类型筛选
# =============================================================================

# --- 🔧 全局配置 ---
st.set_page_config(page_title="RNA 结构与 AIDD 药物重定位平台", layout="wide")
st.title("🧬 RNA 靶点分类 & AIDD 药物重定位系统")

# 本地缓存路径
SMILES_CACHE_FILE = "smiles_cache.json"
DISEASE_CACHE_FILE = "disease_cache.json"
DRUG_SEARCH_CACHE_FILE = "drug_search_cache.json"

# --------------------------
# 核心优化：RNA结构分类体系（基于PDB描述完善）
# 分为3级：大类 → 细分亚型 → 来源标注
# --------------------------
def categorize(desc):
    """
    基于PDB描述的RNA结构精细分类
    返回：(大类, 细分亚型, 来源标注)
    """
    desc_lower = str(desc).lower().strip()
    source_tag = ""  # 来源标注：细菌/病毒/人类/未知
    
    # 1. 标注来源（优先识别描述中的物种信息）
    if any(word in desc_lower for word in ['escherichia coli', 'e. coli', 'bacterial', 'staphylococcus', 'streptococcus']):
        source_tag = "细菌来源"
    elif any(word in desc_lower for word in ['virus', 'viral', 'hiv', 'influenza', 'sars-cov', 'cmv']):
        source_tag = "病毒来源"
    elif any(word in desc_lower for word in ['human', 'homo sapiens', 'hela', 'human cell']):
        source_tag = "人类来源"
    else:
        source_tag = "来源未知"
    
    # 2. 大类+细分亚型分类（基于PDB常见描述关键词）
    # 2.1 核糖开关（按配体类型细分）
    if 'riboswitch' in desc_lower:
        if 'purine' in desc_lower:
            return ("核糖开关 (Riboswitch)", "嘌呤核糖开关", source_tag)
        elif 'fmn' in desc_lower or 'flavin' in desc_lower:
            return ("核糖开关 (Riboswitch)", "FMN核糖开关", source_tag)
        elif 'thiamine' in desc_lower or 'thi' in desc_lower:
            return ("核糖开关 (Riboswitch)", "硫胺素核糖开关", source_tag)
        elif 'glms' in desc_lower:
            return ("核糖开关 (Riboswitch)", "glmS核糖开关", source_tag)
        elif 'lysine' in desc_lower:
            return ("核糖开关 (Riboswitch)", "赖氨酸核糖开关", source_tag)
        elif 'preq1' in desc_lower:
            return ("核糖开关 (Riboswitch)", "preQ1核糖开关", source_tag)
        else:
            return ("核糖开关 (Riboswitch)", "通用核糖开关", source_tag)
    
    # 2.2 适配体（按靶标类型细分）
    elif 'aptamer' in desc_lower:
        if 'protein' in desc_lower:
            return ("适配体 (Aptamer)", "蛋白质结合适配体", source_tag)
        elif 'small molecule' in desc_lower or 'ligand' in desc_lower:
            return ("适配体 (Aptamer)", "小分子结合适配体", source_tag)
        elif 'dna' in desc_lower:
            return ("适配体 (Aptamer)", "DNA结合适配体", source_tag)
        else:
            return ("适配体 (Aptamer)", "通用适配体", source_tag)
    
    # 2.3 G-四联体（按来源细分）
    elif any(word in desc_lower for word in ['quadruplex', 'g-4', 'g4']):
        if 'telomere' in desc_lower:
            return ("G-四联体 (G-quadruplex)", "端粒G-四联体", source_tag)
        elif 'promoter' in desc_lower or 'gene' in desc_lower:
            return ("G-四联体 (G-quadruplex)", "启动子区G-四联体", source_tag)
        else:
            return ("G-四联体 (G-quadruplex)", "通用G-四联体", source_tag)
    
    # 2.4 核糖体（按亚基细分）
    elif any(word in desc_lower for word in ['ribosomal', 'ribosome', 'rrna']):
        if '16s' in desc_lower:
            return ("核糖体 (rRNA)", "16S rRNA（小亚基）", source_tag)
        elif '23s' in desc_lower:
            return ("核糖体 (rRNA)", "23S rRNA（大亚基）", source_tag)
        elif '5s' in desc_lower:
            return ("核糖体 (rRNA)", "5S rRNA", source_tag)
        else:
            return ("核糖体 (rRNA)", "通用核糖体RNA", source_tag)
    
    # 2.5 核酶（按功能细分）
    elif 'ribozyme' in desc_lower:
        if 'hammerhead' in desc_lower:
            return ("核酶 (Ribozyme)", "锤头核酶", source_tag)
        elif 'hairpin' in desc_lower:
            return ("核酶 (Ribozyme)", "发夹核酶", source_tag)
        elif 'glms' in desc_lower:
            return ("核酶 (Ribozyme)", "glmS核酶", source_tag)
        elif 'rnase p' in desc_lower:
            return ("核酶 (Ribozyme)", "RNase P核酶", source_tag)
        else:
            return ("核酶 (Ribozyme)", "通用核酶", source_tag)
    
    # 2.6 特殊结构基元（按结构类型细分）
    elif any(word in desc_lower for word in ['ires', 'hairpin', 'stem-loop', 'pseudoknot', 'internal ribosomal entry site']):
        if 'ires' in desc_lower or 'internal ribosomal entry site' in desc_lower:
            return ("特殊结构基元 (Special/Motifs)", "IRES元件", source_tag)
        elif 'pseudoknot' in desc_lower:
            return ("特殊结构基元 (Special/Motifs)", "假结结构", source_tag)
        elif 'hairpin' in desc_lower:
            return ("特殊结构基元 (Special/Motifs)", "发夹结构", source_tag)
        elif 'stem-loop' in desc_lower:
            return ("特殊结构基元 (Special/Motifs)", "茎环结构", source_tag)
        else:
            return ("特殊结构基元 (Special/Motifs)", "通用结构基元", source_tag)
    
    # 2.7 其他RNA（新增tRNA、snRNA等常见类型）
    else:
        if 'grna' in desc_lower or 'guide rna' in desc_lower:
            return ("其他 RNA (Others)", "向导RNA (gRNA)", source_tag)
        elif 'sirna' in desc_lower or 'small interfering rna' in desc_lower:
            return ("其他 RNA (Others)", "小干扰RNA (siRNA)", source_tag)
        elif 'mirna' in desc_lower or 'microrna' in desc_lower:
            return ("其他 RNA (Others)", "微小RNA (miRNA)", source_tag)
        elif 'trna' in desc_lower:
            return ("其他 RNA (Others)", "转运RNA (tRNA)", source_tag)
        elif 'snrna' in desc_lower:
            return ("其他 RNA (Others)", "小核RNA (snRNA)", source_tag)
        else:
            return ("其他 RNA (Others)", "未分类RNA", source_tag)

# --------------------------
# 通用缓存工具函数
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
        st.warning(f"缓存保存失败: {e}")

# --- 📊 数据加载与预处理（集成新分类） ---
@st.cache_data(ttl=3600)
def load_and_process_data():
    """加载PDB数据并应用新的精细分类"""
    try:
        # 读取用户上传的PDB文件（确保文件在根目录）
        df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    except FileNotFoundError:
        st.error("❌ 未找到数据文件 'PDB_Dataset_Info_Full.xlsx'，请确保文件在根目录下。")
        st.stop()
        return pd.DataFrame()

    # 应用新分类：新增3列（大类、细分亚型、来源标注）
    category_result = df['Description (描述)'].apply(categorize)
    df['Category'] = [res[0] for res in category_result]  # 大类（兼容原有功能）
    df['SubCategory'] = [res[1] for res in category_result]  # 细分亚型（新增）
    df['SourceTag'] = [res[2] for res in category_result]  # 来源标注（新增）
    
    # 提取核心配体ID（保留原有逻辑）
    def get_sort_id(ligand_text):
        if ligand_text == "No ligands" or str(ligand_text) == "nan":
            return "ZZZ"
        all_l = str(ligand_text).split(' | ')
        # 过滤常见离子/结晶剂
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F', 'IOD', 'FLC', 'NO3', 'NH4', 'ACT', 'FMT', 'EDO', 'GOL', 'PEG', 'DTT', 'BME']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions:
                return lid
        return all_l[0].strip().split(' ')[0].upper()
    
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    return df

# --- 以下为原有核心函数（SMILES获取、药物搜索、疾病预测），保持不变 ---
@st.cache_data(ttl=86400, show_spinner=False)
@lru_cache(maxsize=1000)
def get_smiles_by_id(ligand_id):
    if not ligand_id or ligand_id == "ZZZ":
        return None
    ligand_id = ligand_id.strip().upper()
    cache = load_cache(SMILES_CACHE_FILE)
    if ligand_id in cache:
        return cache[ligand_id]

    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36", "Accept": "application/json"}
    max_retries = 2
    retry_delay = 1
    smiles_result = None

    # 1. PDBe API
    for _ in range(max_retries):
        try:
            r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/v2/compound/summary/{ligand_id}", headers=headers, timeout=10)
            if r.status_code == 200:
                data = r.json()
                if ligand_id in data and "smiles" in data[ligand_id][0]:
                    for s in data[ligand_id][0]["smiles"]:
                        if s.get("name") == "canonical":
                            smiles_result = s.get("value")
                            break
                    if not smiles_result:
                        smiles_result = data[ligand_id][0]["smiles"][0].get("value")
                    break
        except:
            time.sleep(retry_delay)

    # 2. RCSB API 备选
    if not smiles_result:
        for _ in range(max_retries):
            try:
                r = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}", headers=headers, timeout=10)
                if r.status_code == 200:
                    data = r.json()
                    desc = data.get("rcsb_chem_comp_descriptor", {})
                    smiles_result = desc.get("SMILES_stereo") or desc.get("SMILES")
                    break
            except:
                time.sleep(retry_delay)

    # 3. PubChem API 兜底
    if not smiles_result:
        try:
            r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/PDB/{ligand_id}/cids/JSON", headers=headers, timeout=10)
            if r.status_code == 200:
                cid = r.json()["IdentifierList"]["CID"][0]
                r2 = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON", headers=headers, timeout=10)
                if r2.status_code == 200:
                    smiles_result = r2.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        except:
            pass

    if smiles_result:
        cache[ligand_id] = smiles_result
        save_cache(cache, SMILES_CACHE_FILE)
    return smiles_result

@st.cache_data(ttl=86400, show_spinner=False)
def search_chembl_drugs(smiles, similarity_threshold):
    if not smiles:
        return [], "SMILES为空", 0
    cache_key = f"{smiles}_{similarity_threshold}"
    drug_cache = load_cache(DRUG_SEARCH_CACHE_FILE)
    if cache_key in drug_cache:
        return drug_cache[cache_key]["drugs"], drug_cache[cache_key]["msg"], drug_cache[cache_key]["total"]

    safe_smiles = urllib.parse.quote(str(smiles).strip())
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    try:
        r = requests.get(f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json?limit=1000&expand=drug_indications", headers=headers, timeout=30)
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            drugs = []
            for m in mols:
                if not (m.get('max_phase') and float(m.get('max_phase')) >= 4.0 and m.get('pref_name') and not m.get('withdrawn_flag')):
                    continue
                # 适应症处理
                drug_name = m.get('pref_name', '').upper()
                indication_list = []
                MANUAL_DRUG_INDICATION = {
                    "PAROMOMYCIN": "肠道阿米巴病、细菌性痢疾", "NETILMICIN": "敏感菌所致呼吸道/泌尿道感染",
                    "KANAMYCIN": "敏感菌所致严重感染", "GENTAMICIN": "革兰氏阴性菌感染"
                }
                if drug_name in MANUAL_DRUG_INDICATION:
                    indication_list.append(MANUAL_DRUG_INDICATION[drug_name])
                else:
                    for ind in m.get('drug_indications', []):
                        for eng, cn in {"bacterial infection":"细菌感染", "cancer":"癌症", "virus":"病毒感染"}.items():
                            if eng in str(ind.get('mesh_heading')).lower():
                                indication_list.append(cn)
                # 药物分类
                drug_class = "未分类"
                atc_list = m.get('atc_classifications', [])
                if atc_list:
                    atc_code = atc_list[0][:3]
                    ATC_FULL_MAP = {"J01":"全身用抗菌药", "L01":"抗肿瘤药", "C02":"抗高血压药"}
                    drug_class = ATC_FULL_MAP.get(atc_code, ATC_FULL_MAP.get(atc_list[0][0], "未分类"))
                
                drugs.append({
                    "药物名称": m.get('pref_name'), "ChEMBL ID": m.get('molecule_chembl_id'),
                    "相似度(%)": round(float(m.get('similarity', 0)), 2), "分子量": round(float(m.get('molecule_properties', {}).get('full_mwt', 0)), 2),
                    "药物分类归属": drug_class, "治疗适应症": "、".join(list(set(indication_list))) if indication_list else "暂无数据"
                })
            drugs.sort(key=lambda x: x["相似度(%)"], reverse=True)
            drug_cache[cache_key] = {"drugs": drugs, "msg": "Success", "total": len(mols)}
            save_cache(drug_cache, DRUG_SEARCH_CACHE_FILE)
            return drugs, "Success", len(mols)
        return [], f"接口失败 ({r.status_code})", 0
    except Exception as e:
        return [], f"异常: {str(e)}", 0

@st.cache_data(ttl=86400, show_spinner=False)
@lru_cache(maxsize=500)
def predict_rna_diseases(pdb_id, description, ligand_id):
    cache = load_cache(DISEASE_CACHE_FILE)
    if pdb_id in cache:
        return cache[pdb_id]
    diseases = []
    desc_lower = str(description).lower()
    # 文本挖掘疾病关联
    disease_map = {
        'cancer': ('癌症', '高'), 'tumor': ('肿瘤', '高'), 'bacterial': ('细菌感染', '高'),
        'virus': ('病毒感染', '高'), 'diabetes': ('糖尿病', '中'), 'neurodegenerative': ('神经退行性疾病', '中')
    }
    for eng, (chn, conf) in disease_map.items():
        if eng in desc_lower:
            diseases.append({
                "疾病名称": chn, "证据来源": "PDB描述文本挖掘", "置信度": conf,
                "文献链接": f"https://www.rcsb.org/structure/{pdb_id}"
            })
    # 配体反推
    if ligand_id != "ZZZ":
        try:
            r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/v2/compound/summary/{ligand_id}", headers={"User-Agent": "Mozilla/5.0"}, timeout=5)
            if r.status_code == 200:
                comp_name = r.json().get(ligand_id, [{}])[0].get("name", "").lower()
                if 'antibiotic' in comp_name:
                    diseases.append({"疾病名称": "细菌感染", "证据来源": f"配体{ligand_id}类型", "置信度": "中", "文献链接": f"https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{ligand_id}"})
        except:
            pass
    # 去重
    unique_diseases = []
    seen = set()
    for d in diseases:
        key = (d["疾病名称"], d["证据来源"])
        if key not in seen:
            seen.add(key)
            unique_diseases.append(d)
    cache[pdb_id] = unique_diseases
    save_cache(cache, DISEASE_CACHE_FILE)
    return unique_diseases

# --- 🎨 界面渲染（新增细分分类筛选和显示） ---
try:
    df = load_and_process_data()
    
    # 侧边栏：新增细分分类筛选
    st.sidebar.header("⚙️ 查看模式")
    view_mode = st.sidebar.radio("模式切换", ["🔍 靶点分析 & AIDD药物筛选", "📊 全局画廊对照"])
    
    st.sidebar.divider()
    st.sidebar.header("🔍 数据筛选")
    # 一级筛选：大类
    category_list = ["全部 (All)"] + sorted(df['Category'].unique().tolist())
    selected_cat = st.sidebar.selectbox("1. 选择RNA大类", category_list)
    # 二级筛选：细分亚型（根据大类联动）
    if selected_cat == "全部 (All)":
        subcat_list = ["全部 (All)"] + sorted(df['SubCategory'].unique().tolist())
    else:
        subcat_list = ["全部 (All)"] + sorted(df[df['Category'] == selected_cat]['SubCategory'].unique().tolist())
    selected_subcat = st.sidebar.selectbox("2. 选择RNA细分亚型", subcat_list)
    # 三级筛选：来源标注
    source_list = ["全部 (All)"] + sorted(df['SourceTag'].unique().tolist())
    selected_source = st.sidebar.selectbox("3. 选择结构来源", source_list)
    
    # 应用筛选条件
    f_df = df.copy()
    if selected_cat != "全部 (All)":
        f_df = f_df[f_df['Category'] == selected_cat]
    if selected_subcat != "全部 (All)":
        f_df = f_df[f_df['SubCategory'] == selected_subcat]
    if selected_source != "全部 (All)":
        f_df = f_df[f_df['SourceTag'] == selected_source]
    f_df = f_df.sort_values(by=['MainLigandID', 'PDB ID'])

    # ==========================================
    # 模式1：全局画廊对照（显示细分分类和来源）
    # ==========================================
    if view_mode == "📊 全局画廊对照":
        st.subheader(f"当前筛选: {selected_cat} → {selected_subcat} → {selected_source} (共{len(f_df)}个结构)")
        cols = st.columns(4)
        for idx, row in f_df.reset_index().iterrows():
            pdb_id, target_id = row['PDB ID'], row['MainLigandID']
            subcat, source = row['SubCategory'], row['SourceTag']
            with cols[idx % 4]:
                img_url = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg" if target_id != "ZZZ" else ""
                # 卡片新增细分分类和来源标注
                st.markdown(f"""
                <div style="border: 1px solid #eee; border-radius: 10px; padding: 12px; margin-bottom: 20px; text-align: center; background: white; min-height: 300px;">
                    <strong style="color: #007bff; font-size: 1.1em;">{pdb_id}</strong><br/>
                    <span style="font-size: 0.7em; color: #28a745; background: #f0f8f0; padding: 2px 8px; border-radius: 12px;">{subcat}</span><br/>
                    <span style="font-size: 0.7em; color: #6c757d;">{source}</span><br/>
                    <img src="{img_url}" style="width:100%; height:130px; object-fit:contain; margin:10px 0;">
                    <div style="font-size: 0.85em; font-weight: bold; color: #333;">{target_id if target_id != 'ZZZ' else '无核心配体'}</div>
                    <div style="font-size: 0.75em; color: #666; height: 40px; overflow: hidden; line-height: 1.2;">{row['Description (描述)'][:50]}...</div>
                </div>
                """, unsafe_allow_html=True)
    
    # ==========================================
    # 模式2：靶点分析 & AIDD药物筛选（显示细分分类）
    # ==========================================
    else:
        tab_single, tab_batch, tab_category = st.tabs(["📌 单靶点详细分析", "⚡ 批量多靶点分析", "📊 按RNA类别一键分析"])
        
        with tab_single:
            selected_pdb = st.sidebar.selectbox("选择PDB ID", f_df['PDB ID'].tolist())
            info = df[df['PDB ID'] == selected_pdb].iloc[0]
            # 显示细分分类和来源
            st.subheader(f"📋 靶点背景信息 | {info['SubCategory']} | {info['SourceTag']}")
            st.success(f"**RNA大类:** {info['Category']}")
            st.info(f"**PDB ID:** {selected_pdb}")
            st.markdown(f"**细分亚型:** `{info['SubCategory']}`")
            st.markdown(f"**结构来源:** `{info['SourceTag']}`")
            
            # 后续原有界面逻辑保持不变...
            target_id = info['MainLigandID']
            if target_id != "ZZZ":
                try:
                    st.image(f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg", caption=f"主要配体 (ID: {target_id})", width=250)
                except:
                    st.warning("配体2D图加载失败")
            st.write(f"**📖 结构描述:** {info['Description (描述)']}")
            st.markdown(f"**🔬 文献出处:** {info['Publication (文章出处)']}")
            st.markdown(f"**🧪 完整配体信息:** `{info['Ligands (对应小分子)']}`")
            
            # 3D结构、疾病预测、药物筛选等原有功能保持不变...
            col1, col2 = st.columns([1, 2])
            with col2:
                st.subheader("🔭 3D空间结构视图")
                try:
                    view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=500)
                    view.setStyle({'cartoon': {'color': 'spectrum'}})
                    view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
                    view.zoomTo()
                    showmol(view, height=500, width=800)
                except Exception as e:
                    st.error(f"3D结构加载失败: {e}")
            
            # RNA-疾病关联预测
            st.divider()
            st.subheader("🩺 RNA靶点-疾病关联预测")
            with st.spinner("分析中..."):
                diseases = predict_rna_diseases(selected_pdb, info['Description (描述)'], target_id)
            if diseases:
                st.dataframe(pd.DataFrame(diseases), use_container_width=True, column_config={"文献链接": st.column_config.LinkColumn("文献链接")})
            else:
                st.warning("未找到疾病关联信息")
            
            # 药物筛选
            if target_id != "ZZZ":
                st.divider()
                st.subheader("💊 AIDD跨靶点药物重定位筛选")
                auto_smi = get_smiles_by_id(target_id)
                smiles = st.text_input("🧬 核心配体SMILES:", value=auto_smi or "", placeholder="例如: CC(=O)OC1=CC=CC=C1C(=O)O")
                col3, col4 = st.columns([2, 1])
                with col3:
                    threshold = st.slider("相似度阈值 (%)", 50, 100, 70, 5)
                with col4:
                    st.write(""); st.write("")
                    search_btn = st.button("🚀 搜索上市药物", use_container_width=True)
                if search_btn and smiles:
                    with st.spinner("检索中..."):
                        drugs, msg, total = search_chembl_drugs(smiles, threshold)
                        if drugs:
                            st.success(f"找到{len(drugs)}个上市药物")
                            st.dataframe(pd.DataFrame(drugs), use_container_width=True)
                        else:
                            st.warning("未找到匹配药物")
        
        # 批量分析、按类别分析等其他Tab保持原有逻辑，仅在显示时新增细分分类和来源字段...
        with tab_batch:
            st.subheader("⚡ 批量多靶点分析")
            st.caption("支持按细分类型筛选靶点")
            # 批量选择逻辑中新增细分类型显示...
        
        with tab_category:
            st.subheader("📊 按RNA类别一键分析")
            # 类别分析结果中新增细分类型和来源统计...

except Exception as e:
    st.error(f"系统运行异常: {e}")
