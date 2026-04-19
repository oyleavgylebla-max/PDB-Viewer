import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol
import requests
import urllib.parse

# 1. 网页全局配置
st.set_page_config(page_title="RNA 结构精细化分类与 AIDD 平台", layout="wide")
st.title("🧬 RNA 靶点分类 & AIDD 药物重定位系统")

@st.cache_data
def load_and_process_data():
    """加载 Excel 并进行分类与配体预处理"""
    df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    
    # --- RNA 精细化分类逻辑 ---
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
    
    # --- 提取核心配体用于排序 ---
    def get_sort_id(ligand_text):
        if ligand_text == "No ligands" or str(ligand_text) == "nan": return "ZZZ"
        all_l = str(ligand_text).split(' | ')
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F', 'IOD', 'FLC', 'NO3', 'NH4', 'ACT', 'FMT', 'EDO', 'GOL', 'PEG', 'DTT', 'BME']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    return df

# --- 🚀 AIDD 核心：三引擎自动获取 SMILES (含 PubChem) ---
@st.cache_data
def get_smiles_from_all_sources(ligand_id):
    """依次尝试 PDBe, RCSB, PubChem 获取 SMILES"""
    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
    }
    
    # 1. 尝试欧洲 PDBe
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{ligand_id}"
        r = requests.get(url, headers=headers, timeout=5)
        if r.status_code == 200:
            data = r.json()
            if ligand_id in data:
                smiles_list = data[ligand_id][0].get("smiles", [])
                for s in smiles_list:
                    if s.get("name") == "canonical": return s.get("value")
                if smiles_list: return smiles_list[0].get("value")
    except: pass

    # 2. 尝试美国 RCSB PDB
    try:
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
        r = requests.get(url, headers=headers, timeout=5)
        if r.status_code == 200:
            data = r.json()
            for desc in data.get("rcsb_chem_comp_descriptor", []):
                if "SMILES" in desc.get("type", "").upper(): return desc.get("descriptor")
    except: pass

    # 3. 🚀 强力后援：PubChem 数据库
    try:
        # 使用配体 ID 作为名称搜索
        pub_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand_id}/property/CanonicalSMILES/JSON"
        r = requests.get(pub_url, headers=headers, timeout=5)
        if r.status_code == 200:
            p_data = r.json()
            smiles = p_data.get("PropertyTable", {}).get("Properties", [{}])[0].get("CanonicalSMILES")
            if smiles: return smiles
    except: pass

    return None

def search_chembl_drugs(smiles, similarity_threshold):
    """在 ChEMBL 中搜索相似 FDA 药物"""
    if not smiles or str(smiles).strip() == "": return [], "SMILES 序列为空"
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json"
    headers = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36", "Accept": "application/json"}
    try:
        r = requests.get(url, headers=headers, timeout=25)
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            drugs = []
            for m in mols:
                phase = m.get('max_phase', 0)
                if phase and float(phase) >= 4.0 and m.get('pref_name'):
                    drugs.append({
                        "药物名称 (Drug)": m.get('pref_name'),
                        "ChEMBL ID": m.get('molecule_chembl_id'),
                        "相似度 (%)": m.get('similarity', 'N/A'),
                        "分子量": m.get('molecule_properties', {}).get('full_mwt', 'N/A')
                    })
            return drugs, "Success"
        else: return [], f"ChEMBL 数据库暂时无法访问 (状态码: {r.status_code})"
    except: return [], "检索超时，请稍后重试"

# --- 主程序逻辑 ---
try:
    df = load_and_process_data()
    st.sidebar.header("⚙️ 显示模式")
    view_mode = st.sidebar.radio("模式切换", ["🔍 单体详细查看 & AIDD分析", "📊 同类全局对比 (画廊)"])
    
    st.sidebar.divider()
    st.sidebar.header("🔍 分类筛选")
    category_list = ["全部 (All)"] + sorted(df['Category'].unique().tolist())
    selected_cat = st.sidebar.selectbox("选择 RNA 类型", category_list)
    
    if selected_cat == "全部 (All)":
        filtered_df = df.sort_values(by=['MainLigandID', 'PDB ID'])
    else:
        filtered_df = df[df['Category'] == selected_cat].sort_values(by=['MainLigandID', 'PDB ID'])

    if view_mode == "📊 同类全局对比 (画廊)":
        st.subheader(f"当前视图: {selected_cat} (已按配体排序)")
        cols = st.columns(4)
        for idx, row in filtered_df.reset_index().iterrows():
            pdb_id, target_id, desc = row['PDB ID'], row['MainLigandID'], row['Description (描述)']
            with cols[idx % 4]:
                img_url = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg" if target_id != "ZZZ" else ""
                st.markdown(f"""
                <div style="border: 1px solid #eee; border-radius: 10px; padding: 12px; margin-bottom: 20px; text-align: center; background: white; min-height: 280px;">
                    <strong style="color: #007bff;">{pdb_id}</strong><br/>
                    <img src="{img_url}" style="width:100%; height:140px; object-fit:contain;">
                    <div style="font-size: 0.8em; font-weight: bold;">{target_id}</div>
                    <div style="font-size: 0.7em; color: #666; height: 45px; overflow: hidden;">{desc}</div>
                </div>
                """, unsafe_allow_html=True)

    else:
        selected_pdb = st.sidebar.selectbox("选择 PDB ID", filtered_df['PDB ID'].tolist())
        info = df[df['PDB ID'] == selected_pdb].iloc[0]
        target_id = info['MainLigandID']
        
        col1, col2 = st.columns([1, 2])
        with col1:
            st.subheader("基本信息")
            st.success(f"分类: {info['Category']}")
            st.info(f"PDB ID: {selected_pdb}")
            if target_id != "ZZZ":
                pdbe_url = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg"
                st.image(pdbe_url, caption=f"配体: {target_id}", width=250)
            with st.expander("详细描述"):
                st.write(info['Description (描述)'])

        with col2:
            st.subheader("3D 空间构象")
            view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=500)
            view.setStyle({'cartoon': {'color': 'spectrum'}}) 
            view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
            view.zoomTo()
            showmol(view, height=500, width=800)

        if target_id != "ZZZ":
            st.divider()
            st.subheader("💊 AIDD 靶点与药物重定位分析")
            
            # 自动抓取 SMILES (三引擎驱动)
            with st.spinner(f"正在多源检索 {target_id} 的化学编码..."):
                auto_smiles = get_smiles_from_all_sources(target_id)
            
            if not auto_smiles:
                st.warning(f"⚠️ 自动抓取失败。PDB 与 PubChem 库中均未找到 {target_id} 的规范 SMILES。")
            
            smiles = st.text_input("🧬 配体 SMILES (支持自动抓取结果/手动覆盖):", value=auto_smiles if auto_smiles else "")
            
            if smiles:
                c1, c2 = st.columns([2, 1])
                with c1:
                    sim_threshold = st.slider("Tanimoto 相似度阈值 (%)", 70, 100, 70, 5)
                with c2:
                    st.write("")
                    st.write("")
                    if st.button("🚀 开始跨靶点搜索", use_container_width=True):
                        with st.spinner("正在扫描 ChEMBL 上市药物库..."):
                            fda_drugs, msg = search_chembl_drugs(smiles, sim_threshold)
                            if fda_drugs:
                                st.success(f"🎉 找到 {len(fda_drugs)} 个相似的 FDA 批准药物！")
                                st.dataframe(pd.DataFrame(fda_drugs), use_container_width=True)
                            elif msg == "Success":
                                st.warning("🔍 检索成功，但未匹配到相似度足够高的上市药物。")
                            else:
                                st.error(f"🚨 错误: {msg}")

except Exception as e:
    st.error(f"运行中发生错误: {e}")
