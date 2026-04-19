import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol
import requests
import urllib.parse

st.set_page_config(page_title="RNA 结构精细化分类与 AIDD 平台", layout="wide")
st.title("🧬 RNA 靶点分类 & AIDD 药物重定位系统")

@st.cache_data
def load_and_process_data():
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
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F', 'IOD', 'FLC', 'NO3', 'NH4', 'ACT', 'FMT', 'EDO', 'GOL', 'PEG', 'DTT', 'BME']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    return df

@st.cache_data
def get_smiles_from_pdb_ids(ligand_id):
    headers = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/120.0.0.0 Safari/537.36"}
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
    try:
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
        r = requests.get(url, headers=headers, timeout=5)
        if r.status_code == 200:
            data = r.json()
            for desc in data.get("rcsb_chem_comp_descriptor", []):
                if "SMILES" in desc.get("type", "").upper(): return desc.get("descriptor")
    except: pass
    return None

@st.cache_data
def get_smiles_from_pubchem_name(name):
    if not name or str(name).strip() == "": return None
    clean_name = str(name).strip()
    if clean_name.startswith('[') and clean_name.endswith(']'):
        clean_name = clean_name[1:-1].strip()
    safe_name = urllib.parse.quote(clean_name)
    pub_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{safe_name}/property/CanonicalSMILES/JSON"
    headers = {"User-Agent": "Mozilla/5.0"}
    try:
        r = requests.get(pub_url, headers=headers, timeout=10)
        if r.status_code == 200:
            return r.json().get("PropertyTable", {}).get("Properties", [{}])[0].get("CanonicalSMILES")
    except: pass
    return None

def search_chembl_drugs(smiles, similarity_threshold):
    if not smiles or str(smiles).strip() == "": return [], "SMILES为空", 0
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    # 💡 核心升级：强行让 API 返回前 1000 个相似分子，打破分页陷阱！
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json?limit=1000"
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    try:
        r = requests.get(url, headers=headers, timeout=25)
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            total_found = len(mols) # 记录总共找到了多少个相似分子
            drugs = []
            for m in mols:
                phase = m.get('max_phase', 0)
                # 过滤出 Phase 4 (已上市) 的分子
                if phase and float(phase) >= 4.0 and m.get('pref_name'):
                    drugs.append({
                        "药物名称 (Drug)": m.get('pref_name'),
                        "ChEMBL ID": m.get('molecule_chembl_id'),
                        "相似度 (%)": m.get('similarity', 'N/A'),
                        "分子量": m.get('molecule_properties', {}).get('full_mwt', 'N/A')
                    })
            return drugs, "Success", total_found
        else: return [], f"接口被拒 (状态码: {r.status_code})", 0
    except: return [], "检索超时", 0

# --- 主程序逻辑 ---
try:
    df = load_and_process_data()
    st.sidebar.header("⚙️ 显示模式")
    view_mode = st.sidebar.radio("模式切换", ["🔍 单体详细查看 & AIDD分析", "📊 同类全局对比 (画廊)"])
    st.sidebar.divider()
    category_list = ["全部 (All)"] + sorted(df['Category'].unique().tolist())
    selected_cat = st.sidebar.selectbox("选择 RNA 类型", category_list)
    
    if selected_cat == "全部 (All)": filtered_df = df.sort_values(by=['MainLigandID', 'PDB ID'])
    else: filtered_df = df[df['Category'] == selected_cat].sort_values(by=['MainLigandID', 'PDB ID'])

    if view_mode == "📊 同类全局对比 (画廊)":
        st.subheader(f"当前视图: {selected_cat}")
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
            with st.expander("详细描述"): st.write(info['Description (描述)'])

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
            
            auto_smiles = get_smiles_from_pdb_ids(target_id)
            raw_ligand_str = str(info['Ligands (对应小分子)'])
            name_guess = raw_ligand_str.replace(target_id, "").replace("|", "").replace(";", "").strip()
            if name_guess.startswith('[') and name_guess.endswith(']'): name_guess = name_guess[1:-1].strip()
            
            st.write("### 🪪 数据库全名搜寻通道")
            c1, c2 = st.columns([3, 1])
            with c1: input_name = st.text_input("小分子全名 (自动去除干扰符号):", value=name_guess)
            with c2:
                st.write(""); st.write("")
                if st.button("🌐 去 PubChem 搜全名", use_container_width=True):
                    with st.spinner(f"正在检索 '{input_name}'..."):
                        pub_smiles = get_smiles_from_pubchem_name(input_name)
                        if pub_smiles:
                            auto_smiles = pub_smiles
                            st.success("🎉 检索成功！(注意：过于复杂的长串 IUPAC 名称仍建议使用网页版手动查SMILES)")
                        else:
                            st.error("❌ API 匹配失败。建议直接在 PubChem 网页查出 SMILES 后粘贴到下方。")

            st.write("---")
            smiles = st.text_input("🧬 最终用于靶点匹配的 SMILES 序列:", value=auto_smiles if auto_smiles else "")
            
            if smiles:
                c3, c4 = st.columns([2, 1])
                with c3: sim_threshold = st.slider("Tanimoto 相似度阈值 (%)", 50, 100, 70, 5)
                with c4:
                    st.write(""); st.write("")
                    if st.button("🚀 开始跨靶点搜索 (ChEMBL)", use_container_width=True):
                        with st.spinner("正在提取最多前 1000 个相似分子进行 Phase 4 筛选..."):
                            fda_drugs, msg, total_found = search_chembl_drugs(smiles, sim_threshold)
                            if fda_drugs:
                                st.success(f"🎉 在 {total_found} 个相似化合物中，成功筛选出 {len(fda_drugs)} 个 FDA 批准药物！")
                                st.dataframe(pd.DataFrame(fda_drugs), use_container_width=True)
                            elif msg == "Success":
                                # 💡 新增透视提示：告诉你不是代码没搜到，而是确实没上市药！
                                st.warning(f"🔍 检索完成！在 ChEMBL 中找到了 **{total_found}** 个与它相似度 >{sim_threshold}% 的分子。")
                                st.info("🧪 但是，这其中 **没有任何一个是 FDA 已上市的药物 (Phase 4)**。这些相似分子大多处于实验室合成阶段。")
                            else:
                                st.error(f"🚨 错误: {msg}")

except Exception as e:
    st.error(f"运行中发生错误: {e}")
