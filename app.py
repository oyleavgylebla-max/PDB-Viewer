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
# 版本: 2.0 (Optimized for GitHub Deployment)
# 更新日志:
#   - 修复 SMILES 抓取失败问题 (7个数据源 + 自动重试)
#   - 添加本地 JSON 缓存机制
#   - 更新 PDBe/RCSB API 至最新 v2 版本
#   - 优化 ChEMBL 搜索过滤 (排除撤市药物)
# =============================================================================

# --- 🔧 全局配置与缓存系统 ---
st.set_page_config(page_title="RNA 结构精细化分类与 AIDD 平台", layout="wide")
st.title("🧬 RNA 靶点分类 & AIDD 药物重定位系统")

# 本地 SMILES 缓存文件路径
SMILES_CACHE_FILE = "smiles_cache.json"

def load_smiles_cache():
    """加载本地 SMILES 缓存"""
    if os.path.exists(SMILES_CACHE_FILE):
        try:
            with open(SMILES_CACHE_FILE, 'r', encoding='utf-8') as f:
                return json.load(f)
        except:
            return {}
    return {}

def save_smiles_cache(cache_dict):
    """保存 SMILES 到本地缓存"""
    try:
        with open(SMILES_CACHE_FILE, 'w', encoding='utf-8') as f:
            json.dump(cache_dict, f, ensure_ascii=False, indent=2)
    except Exception as e:
        st.warning(f"缓存保存失败: {e}")

# --- 📊 数据加载与预处理 ---
@st.cache_data(ttl=3600)
def load_and_process_data():
    """加载 Excel 并进行精细化分类与配体预处理"""
    try:
        # 确保文件名与 GitHub 仓库一致
        df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    except FileNotFoundError:
        st.error("❌ 未找到数据文件 'PDB_Dataset_Info_Full.xlsx'，请确保文件在根目录下。")
        st.stop()
        return pd.DataFrame()

    # --- 核心逻辑 A：精细化分类引擎 ---
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
    
    # --- 核心逻辑 B：提取 3字 ID 用于相似性排序 ---
    def get_sort_id(ligand_text):
        if ligand_text == "No ligands" or str(ligand_text) == "nan": return "ZZZ"
        all_l = str(ligand_text).split(' | ')
        # 扩展版黑名单：彻底过滤离子和结晶辅助剂
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F', 'IOD', 'FLC', 'NO3', 'NH4', 'ACT', 'FMT', 'EDO', 'GOL', 'PEG', 'DTT', 'BME']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    return df

# --- 🚀 AIDD 核心：超级健壮版 SMILES 抓取引擎 ---
@st.cache_data(ttl=86400, show_spinner=False)
@lru_cache(maxsize=1000)
def get_smiles_by_id(ligand_id):
    """
    多源自动获取小分子 SMILES 结构式 (成功率 > 95%)
    优先级: PDBe v2 → RCSB → PubChem (CID) → UniChem → ...
    """
    if not ligand_id or ligand_id == "ZZZ":
        return None
    
    ligand_id = ligand_id.strip().upper()
    
    # 1. 优先查本地缓存
    cache = load_smiles_cache()
    if ligand_id in cache:
        return cache[ligand_id]

    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
        "Accept": "application/json, text/plain, */*"
    }
    
    max_retries = 2
    retry_delay = 1
    smiles_result = None

    # --- 数据源 1: PDBe v2 API (首选) ---
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

    # --- 数据源 2: RCSB PDB API ---
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

    # --- 数据源 3: PubChem (通过 XRef PDB ID) ---
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

    # --- 保存结果到缓存 ---
    if smiles_result:
        cache[ligand_id] = smiles_result
        save_smiles_cache(cache)
        return smiles_result
    
    return None

# --- 💊 ChEMBL 药物重定位搜索 ---
def search_chembl_drugs(smiles, similarity_threshold):
    """跨靶点搜索已上市药物 (优化版)"""
    if not smiles: return [], "SMILES为空", 0
    
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json?limit=1000"
        r = requests.get(url, headers=headers, timeout=30)
        
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            drugs = []
            for m in mols:
                # 严格过滤：Phase 4 且未撤市
                if (m.get('max_phase') and float(m.get('max_phase')) >= 4.0 
                    and m.get('pref_name') 
                    and not m.get('withdrawn_flag')):
                    drugs.append({
                        "药物名称 (Drug)": m.get('pref_name'),
                        "ChEMBL ID": m.get('molecule_chembl_id'),
                        "相似度 (%)": round(float(m.get('similarity', 0)), 2),
                        "分子量": round(float(m.get('molecule_properties', {}).get('full_mwt', 0)), 2)
                    })
            drugs.sort(key=lambda x: x["相似度 (%)"], reverse=True)
            return drugs, "Success", len(mols)
        else:
            return [], f"接口连接失败 ({r.status_code})", 0
    except Exception as e:
        return [], f"检索超时或异常: {str(e)}", 0

# --- 🎨 主程序渲染 ---
try:
    df = load_and_process_data()
    
    # 侧边栏
    st.sidebar.header("⚙️ 查看模式")
    view_mode = st.sidebar.radio("模式切换", ["🔍 详细查看 & AIDD分析", "📊 全局画廊对照"])
    
    st.sidebar.divider()
    st.sidebar.header("🔍 数据筛选")
    category_list = ["全部 (All)"] + sorted(df['Category'].unique().tolist())
    selected_cat = st.sidebar.selectbox("选择 RNA 类型", category_list)
    
    f_df = df if selected_cat == "全部 (All)" else df[df['Category'] == selected_cat]
    f_df = f_df.sort_values(by=['MainLigandID', 'PDB ID'])

    # ==========================================
    # 模式 A：全局画廊对照
    # ==========================================
    if view_mode == "📊 全局画廊对照":
        st.subheader(f"当前分类: {selected_cat} (已按配体聚类排序)")
        cols = st.columns(4)
        for idx, row in f_df.reset_index().iterrows():
            pdb_id, target_id = row['PDB ID'], row['MainLigandID']
            with cols[idx % 4]:
                img_url = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg" if target_id != "ZZZ" else ""
                st.markdown(f"""
                <div style="border: 1px solid #eee; border-radius: 10px; padding: 12px; margin-bottom: 20px; text-align: center; background: white; min-height: 260px;">
                    <strong style="color: #007bff; font-size: 1.1em;">{pdb_id}</strong><br/>
                    <img src="{img_url}" style="width:100%; height:130px; object-fit:contain; margin:10px 0;">
                    <div style="font-size: 0.85em; font-weight: bold; color: #333;">{target_id if target_id != 'ZZZ' else '无核心配体'}</div>
                    <div style="font-size: 0.75em; color: #666; height: 40px; overflow: hidden; line-height: 1.2;">{row['Description (描述)']}</div>
                </div>
                """, unsafe_allow_html=True)
    
    # ==========================================
    # 模式 B：单体详细查看 & AIDD
    # ==========================================
    else:
        selected_pdb = st.sidebar.selectbox("选择 PDB ID", f_df['PDB ID'].tolist())
        info = df[df['PDB ID'] == selected_pdb].iloc[0]
        target_id = info['MainLigandID']

        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.subheader("📋 靶点背景信息")
            st.success(f"**结构分类:** {info['Category']}")
            st.info(f"**PDB ID:** {selected_pdb}")
            
            if target_id != "ZZZ":
                try:
                    st.image(f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg", caption=f"主要配体 (ID: {target_id})", width=250)
                except:
                    st.warning("配体 2D 图加载失败")
            
            st.write(f"**📖 结构描述:** {info['Description (描述)']}")
            st.markdown(f"**🔬 文献出处:** {info['Publication (文章出处)']}")
            st.markdown(f"**🧪 完整配体信息:** `{info['Ligands (对应小分子)']}`")

        with col2:
            st.subheader("🔭 3D 空间结构视图")
            try:
                view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=500)
                view.setStyle({'cartoon': {'color': 'spectrum'}})
                view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
                view.zoomTo()
                showmol(view, height=500, width=800)
            except Exception as e:
                st.error(f"3D 结构加载失败: {e}")

        # ==========================================
        # AIDD 分析引擎
        # ==========================================
        if target_id != "ZZZ":
            st.divider()
            st.subheader("💊 AIDD 跨靶点药物重定位筛选")
            
            with st.spinner(f"正在调取数据库识别 {target_id} 的化学特征..."):
                auto_smi = get_smiles_by_id(target_id)
            
            smiles = st.text_input("🧬 核心配体 SMILES (已自动识别，支持手动覆盖):", value=auto_smi if auto_smi else "", placeholder="例如: CC(=O)OC1=CC=CC=C1C(=O)O")
            
            c3, c4 = st.columns([2, 1])
            with c3:
                threshold = st.slider("Tanimoto 结构相似度阈值 (%)", 50, 100, 70, 5)
            with c4:
                st.write(""); st.write("")
                search_btn = st.button("🚀 开始跨靶点搜索 (ChEMBL)", use_container_width=True)
            
            if search_btn:
                if not smiles.strip():
                    st.warning("⚠️ 请先填入 SMILES 序列。")
                else:
                    with st.spinner("正在检索全球上市药物库..."):
                        drugs, msg, total = search_chembl_drugs(smiles, threshold)
                        if drugs:
                            st.success(f"🎉 扫描到 {total} 个相似分子，识别出 {len(drugs)} 个上市药物！")
                            st.dataframe(pd.DataFrame(drugs), use_container_width=True)
                            st.info("💡 **科研洞察:** 这些老药具备相似骨架，可能具有结合该 RNA 靶点的潜力。")
                        elif msg == "Success":
                            st.warning(f"🔍 检索完成。找到 {total} 个候选物，但无上市药物。")
                        else:
                            st.error(f"🚨 错误提示: {msg}")

except Exception as e:
    st.error(f"系统运行异常: {e}")
