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
    """加载数据并进行分类预处理"""
    df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    
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
    
    def get_sort_id(ligand_text):
        if ligand_text == "No ligands" or str(ligand_text) == "nan": return "ZZZ"
        all_l = str(ligand_text).split(' | ')
        # 常见杂质/离子黑名单
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F', 'IOD', 'FLC', 'NO3', 'NH4', 'ACT', 'FMT', 'EDO', 'GOL', 'PEG', 'DTT', 'BME']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    return df

# --- 🚀 AIDD 核心：三引擎自动获取 SMILES (基于 3字 ID) ---
@st.cache_data
def get_smiles_by_id(ligand_id):
    """基于 3字 ID 依次尝试 PDBe, RCSB, PubChem"""
    if not ligand_id or ligand_id == "ZZZ": return None
    
    headers = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/120.0.0.0 Safari/537.36"}
    
    # 1. PDBe 引擎
    try:
        r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{ligand_id}", headers=headers, timeout=5)
        if r.status_code == 200:
            smiles_list = r.json().get(ligand_id, [{}])[0].get("smiles", [])
            for s in smiles_list:
                if s.get("name") == "canonical": return s.get("value")
            if smiles_list: return smiles_list[0].get("value")
    except: pass

    # 2. PubChem 引擎 (使用 ID 匹配名)
    try:
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand_id}/property/CanonicalSMILES/JSON", headers=headers, timeout=5)
        if r.status_code == 200:
            return r.json().get("PropertyTable", {}).get("Properties", [{}])[0].get("CanonicalSMILES")
    except: pass
    
    return None

def search_chembl_drugs(smiles, similarity_threshold):
    """在 ChEMBL 中搜索相似药物 (含 1000 结果突破)"""
    if not smiles: return [], "SMILES为空", 0
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    # 💡 limit=1000 确保不会漏掉排在后面的上市药物
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json?limit=1000"
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    try:
        r = requests.get(url, headers=headers, timeout=25)
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            drugs = []
            for m in mols:
                if m.get('max_phase') and float(m.get('max_phase')) >= 4.0 and m.get('pref_name'):
                    drugs.append({
                        "药物名称 (Drug)": m.get('pref_name'),
                        "ChEMBL ID": m.get('molecule_chembl_id'),
                        "相似度 (%)": m.get('similarity', 'N/A'),
                        "分子量": m.get('molecule_properties', {}).get('full_mwt', 'N/A')
                    })
            return drugs, "Success", len(mols)
        else: return [], f"接口被拒 ({r.status_code})", 0
    except: return [], "检索超时", 0

# --- 主程序渲染 ---
try:
    df = load_and_process_data()
    st.sidebar.header("⚙️ 查看选项")
    view_mode = st.sidebar.radio("模式", ["🔍 详细查看 & AIDD", "📊 全局画廊"])
    
    selected_cat = st.sidebar.selectbox("RNA 分类", ["全部"] + sorted(df['Category'].unique().tolist()))
    f_df = df if selected_cat == "全部" else df[df['Category'] == selected_cat]
    f_df = f_df.sort_values(by=['MainLigandID', 'PDB ID'])

    if view_mode == "📊 全局画廊":
        st.subheader(f"画廊模式 ({selected_cat})")
        cols = st.columns(4)
        for idx, row in f_df.reset_index().iterrows():
            pdb, lid = row['PDB ID'], row['MainLigandID']
            with cols[idx % 4]:
                img = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{lid}_400.svg" if lid != "ZZZ" else ""
                st.markdown(f'<div style="border:1px solid #eee;border-radius:10px;padding:10px;text-align:center;background:white;min-height:250px;"><strong>{pdb}</strong><br/><img src="{img}" style="height:130px;object-fit:contain;"><br/>{lid}</div>', unsafe_allow_html=True)
    
    else:
        selected_pdb = st.sidebar.selectbox("选择 PDB ID", f_df['PDB ID'].tolist())
        info = df[df['PDB ID'] == selected_pdb].iloc[0]
        target_id = info['MainLigandID']

        col1, col2 = st.columns([1, 2])
        with col1:
            st.success(f"分类: {info['Category']}")
            if target_id != "ZZZ":
                st.image(f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg", caption=f"ID: {target_id}", width=250)
            st.write(f"描述: {info['Description (描述)']}")

        with col2:
            view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=500)
            view.setStyle({'cartoon': {'color': 'spectrum'}})
            view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
            view.zoomTo(); showmol(view, height=500, width=800)

        if target_id != "ZZZ":
            st.divider()
            st.subheader("💊 AIDD 跨靶点药物筛选")
            # 自动基于 ID 抓取
            with st.spinner(f"正在抓取 {target_id} 的结构信息..."):
                auto_smi = get_smiles_by_id(target_id)
            
            smiles = st.text_input("🧬 核心配体 SMILES:", value=auto_smi if auto_smi else "")
            
            if smiles:
                c3, c4 = st.columns([2, 1])
                with c3: threshold = st.slider("结构相似度阈值 (%)", 50, 100, 70, 5)
                with c4:
                    st.write(""); st.write("")
                    if st.button("🚀 搜索相似上市药", use_container_width=True):
                        with st.spinner("检索中..."):
                            drugs, msg, total = search_chembl_drugs(smiles, threshold)
                            if drugs:
                                st.success(f"🎉 在 {total} 个相似分子中揪出 {len(drugs)} 个上市药物！")
                                st.dataframe(pd.DataFrame(drugs), use_container_width=True)
                            elif msg == "Success":
                                st.warning(f"🔍 检索完成。找到 {total} 个相似分子，但其中没有 FDA 已上市药物。")
                            else: st.error(f"🚨 错误: {msg}")

except Exception as e:
    st.error(f"发生错误: {e}")
