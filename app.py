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
    """加载 Excel 并进行精细化分类与配体预处理"""
    # 确保文件名与 GitHub 仓库一致
    df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    
    # --- 核心逻辑 A：精细化分类引擎 (含 G4 和 rRNA) ---
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

# --- 🚀 AIDD 核心：三引擎自动获取 SMILES (基于 3字 ID) ---
@st.cache_data
def get_smiles_by_id(ligand_id):
    """依次尝试 PDBe, RCSB, PubChem 抓取，带浏览器伪装"""
    if not ligand_id or ligand_id == "ZZZ": return None
    headers = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/120.0.0.0 Safari/537.36"}
    
    # 1. 尝试 PDBe (最稳)
    try:
        r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{ligand_id}", headers=headers, timeout=5)
        if r.status_code == 200:
            smiles_list = r.json().get(ligand_id, [{}])[0].get("smiles", [])
            for s in smiles_list:
                if s.get("name") == "canonical": return s.get("value")
            if smiles_list: return smiles_list[0].get("value")
    except: pass

    # 2. 尝试美国 RCSB
    try:
        r = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}", headers=headers, timeout=5)
        if r.status_code == 200:
            for desc in r.json().get("rcsb_chem_comp_descriptor", []):
                if "SMILES" in desc.get("type", "").upper(): return desc.get("descriptor")
    except: pass

    # 3. 尝试 PubChem (后援)
    try:
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand_id}/property/CanonicalSMILES/JSON", headers=headers, timeout=5)
        if r.status_code == 200:
            return r.json().get("PropertyTable", {}).get("Properties", [{}])[0].get("CanonicalSMILES")
    except: pass
    return None

def search_chembl_drugs(smiles, similarity_threshold):
    """跨靶点搜索已上市药物 (带 1000 条深度穿透)"""
    if not smiles: return [], "SMILES为空", 0
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    # 💡 穿透分页陷阱：limit=1000
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json?limit=1000"
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    try:
        r = requests.get(url, headers=headers, timeout=25)
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            drugs = []
            for m in mols:
                # 过滤已上市药物 (Phase 4)
                if m.get('max_phase') and float(m.get('max_phase')) >= 4.0 and m.get('pref_name'):
                    drugs.append({
                        "药物名称 (Drug)": m.get('pref_name'),
                        "ChEMBL ID": m.get('molecule_chembl_id'),
                        "相似度 (%)": m.get('similarity', 'N/A'),
                        "分子量": m.get('molecule_properties', {}).get('full_mwt', 'N/A')
                    })
            return drugs, "Success", len(mols)
        else: return [], f"接口连接失败 ({r.status_code})", 0
    except: return [], "检索超时", 0

# --- 主程序渲染 ---
try:
    df = load_and_process_data()
    
    # 2. 侧边栏
    st.sidebar.header("⚙️ 查看模式")
    view_mode = st.sidebar.radio("模式切换", ["🔍 详细查看 & AIDD分析", "📊 全局画廊对照"])
    
    st.sidebar.divider()
    st.sidebar.header("🔍 数据筛选")
    category_list = ["全部 (All)"] + sorted(df['Category'].unique().tolist())
    selected_cat = st.sidebar.selectbox("选择 RNA 类型", category_list)
    
    # 筛选并按照配体 ID 排序（实现相似分子聚类显示）
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
    # 模式 B：单体详细查看 (重点修复补全区)
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
            
            # 展示核心配体 2D 结构
            if target_id != "ZZZ":
                st.image(f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg", caption=f"主要配体 (ID: {target_id})", width=250)
            
            # --- 💡 重点修复：补全缺失字段 ---
            st.write(f"**📖 结构描述:** {info['Description (描述)']}")
            st.markdown(f"**🔬 文献出处:** {info['Publication (文章出处)']}")
            st.markdown(f"**🧪 完整配体信息:** `{info['Ligands (对应小分子)']}`")

        with col2:
            st.subheader("🔭 3D 空间结构视图")
            view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=500)
            view.setStyle({'cartoon': {'color': 'spectrum'}})
            view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
            view.zoomTo()
            showmol(view, height=500, width=800)

        # ==========================================
        # 🚀 AIDD 分析引擎 (UI 逻辑修复区)
        # ==========================================
        if target_id != "ZZZ":
            st.divider()
            st.subheader("💊 AIDD 跨靶点药物重定位筛选")
            
            # 自动抓取 SMILES
            with st.spinner(f"正在调取数据库识别 {target_id} 的化学特征..."):
                auto_smi = get_smiles_by_id(target_id)
            
            # 输入框（即使没抓到也保留，方便手动粘贴）
            smiles = st.text_input("🧬 核心配体 SMILES (已自动识别，支持手动覆盖):", value=auto_smi if auto_smi else "")
            
            # 💡 核心修复：滑动条和按钮改为全时显示，不受 smiles 变量值限制
            c3, c4 = st.columns([2, 1])
            with c3:
                threshold = st.slider("Tanimoto 结构相似度阈值 (%)", 50, 100, 70, 5, help="推荐 70-80% 以获得高活性类似物")
            with c4:
                st.write(""); st.write("") # 对齐占位
                search_btn = st.button("🚀 开始跨靶点搜索 (ChEMBL)", use_container_width=True)
            
            # 执行搜索逻辑
            if search_btn:
                if not smiles.strip():
                    st.warning("⚠️ 请先在输入框中填入分子的 SMILES 序列（例如阿司匹林: CC(=O)OC1=CC=CC=C1C(=O)O）。")
                else:
                    with st.spinner("正在检索全球上市药物库 (扫描 Phase 4 候选物)..."):
                        drugs, msg, total = search_chembl_drugs(smiles, threshold)
                        if drugs:
                            st.success(f"🎉 在扫描到的 {total} 个相似分子中，精准识别出 {len(drugs)} 个 FDA 上市药物！")
                            st.dataframe(pd.DataFrame(drugs), use_container_width=True)
                            st.info("💡 **科研洞察:** 这些老药具备相似的化学骨架，可能具有结合该 RNA 靶点的潜力。")
                        elif msg == "Success":
                            st.warning(f"🔍 检索完成。找到 {total} 个相似候选物，但无一属于 FDA 上市药物 (Phase 4)。")
                        else:
                            st.error(f"🚨 错误提示: {msg}")

except Exception as e:
    st.error(f"系统运行中发生异常: {e}")
