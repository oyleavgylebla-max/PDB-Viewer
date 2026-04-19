import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol
import requests
import urllib.parse

# 1. 网页配置
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
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'BR', 'I', 'F', 'EDO', 'GOL']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    return df

# --- 新增 AIDD 功能函数 ---
@st.cache_data
def get_smiles_from_pdb(ligand_id):
    """从 PDB 数据库获取分子的 SMILES 表达式"""
    try:
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
        r = requests.get(url, timeout=5)
        if r.status_code == 200:
            data = r.json()
            for desc in data.get("rcsb_chem_comp_descriptor", []):
                if desc.get("type") == "SMILES":
                    return desc.get("descriptor")
    except:
        return None
    return None

def search_chembl_drugs(smiles, similarity_threshold):
    """在 ChEMBL 中搜索相似分子，并过滤出 FDA 批准药物 (Max Phase = 4)"""
    safe_smiles = urllib.parse.quote(smiles)
    # ChEMBL 相似度 API
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json"
    try:
        r = requests.get(url, headers={"Accept": "application/json"}, timeout=15)
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            # 过滤：只要临床 4 期（已上市药物），且有名字的
            drugs = []
            for m in mols:
                if m.get('max_phase') == 4 and m.get('pref_name'):
                    drugs.append({
                        "药物名称 (Drug)": m.get('pref_name'),
                        "ChEMBL ID": m.get('molecule_chembl_id'),
                        "相似度 (%)": m.get('similarity', 'N/A'),
                        "分子式": m.get('molecule_properties', {}).get('full_mwt', 'N/A')
                    })
            return drugs
    except:
        return []
    return []

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

    # ==========================================
    # 模式 B：同类全局对比 (画廊模式) - 保持不变
    # ==========================================
    if view_mode == "📊 同类全局对比 (画廊)":
        st.subheader(f"当前视图: {selected_cat} (已按配体相似性排序)")
        cols_per_row = 4
        columns = st.columns(cols_per_row)
        for idx, row in filtered_df.reset_index().iterrows():
            col_idx = idx % cols_per_row
            pdb_id = row['PDB ID']
            target_id = row['MainLigandID']
            desc = row['Description (描述)']
            with columns[col_idx]:
                if target_id != "ZZZ":
                    img_url = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg"
                    img_html = f'<img src="{img_url}" style="width:100%; height:140px; object-fit:contain; margin-bottom:5px;">'
                else:
                    img_html = '<div style="height:140px; display:flex; align-items:center; justify-content:center; color:#ccc; background:#f9f9f9; border-radius:5px; margin-bottom:5px;">无核心配体</div>'
                card_html = f"""
                <div style="border: 1px solid #eee; border-radius: 10px; padding: 12px; margin-bottom: 20px; text-align: center; background: white;">
                    <strong style="color: #007bff; font-size: 1.1em;">{pdb_id}</strong>
                    {img_html}
                    <div style="font-size: 0.85em; color: #333; font-weight: bold; margin-bottom: 5px;">{target_id if target_id != 'ZZZ' else '-'}</div>
                    <div style="font-size: 0.75em; color: #666; height: 45px; overflow: hidden;">{desc}</div>
                </div>
                """
                st.markdown(card_html, unsafe_allow_html=True)

    # ==========================================
    # 模式 A：单体详细查看 & AIDD 靶点发现
    # ==========================================
    else:
        selected_pdb = st.sidebar.selectbox("选择 PDB ID", filtered_df['PDB ID'].tolist())
        info = df[df['PDB ID'] == selected_pdb].iloc[0]
        target_id = info['MainLigandID']
        
        col1, col2 = st.columns([1, 2])
        with col1:
            st.subheader("基本信息")
            st.success(f"**分类:** {info['Category']}")
            st.info(f"**PDB ID:** {selected_pdb}")
            if target_id != "ZZZ":
                st.markdown(f"**核心配体:** `{target_id}`")
                pdbe_url = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg"
                st.image(pdbe_url, width=200)
            
            with st.expander("查看完整描述"):
                st.write(info['Description (描述)'])
                st.write(f"文献: {info['Publication (文章出处)']}")

        with col2:
            st.subheader("3D 空间构象")
            view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=450)
            view.setStyle({'cartoon': {'color': 'spectrum'}}) 
            view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
            view.zoomTo()
            showmol(view, height=450, width=800)

        # -----------------------------------------
        # 🚀 新增模块：AIDD 药物重定位分析
        # -----------------------------------------
        if target_id != "ZZZ":
            st.divider()
            st.subheader("💊 AIDD 靶点与药物重定位分析 (Drug Repurposing)")
            st.write("利用 ChEMBL 数据库，寻找与当前 RNA 天然配体结构相似的 **FDA 已批准上市药物**。")
            
            # 1. 提取 SMILES
            smiles = get_smiles_from_pdb(target_id)
            if smiles:
                st.code(f"天然配体 SMILES: {smiles}", language="text")
                
                # 2. 设置相似度阈值和搜索按钮
                c1, c2 = st.columns([2, 1])
                with c1:
                    sim_threshold = st.slider("Tanimoto 结构相似度阈值 (%)", min_value=50, max_value=100, value=70, step=5, 
                                            help=">80% 通常具有高度相似活性。调低至 50% 可以发现更多潜在的跨靶点母核结构。")
                with c2:
                    st.write("") # 占位对齐
                    st.write("")
                    search_btn = st.button("🚀 搜索相似 FDA 药物", use_container_width=True)
                
                # 3. 执行搜索
                if search_btn:
                    with st.spinner(f"正在扫描 ChEMBL 数据库 (寻找相似度 > {sim_threshold}% 的上市药物)..."):
                        fda_drugs = search_chembl_drugs(smiles, sim_threshold)
                        
                        if fda_drugs:
                            st.success(f"🎉 发现 {len(fda_drugs)} 个相似的 FDA 批准药物！")
                            st.dataframe(pd.DataFrame(fda_drugs), use_container_width=True)
                            st.info("💡 **科研洞察:** 如果以上药物原本是治疗某疾病的（如抗癌、抗菌），说明该小分子可能通过脱靶效应结合了当前的 RNA 结构；或者，这个 RNA 本身就是该疾病的潜在治疗靶点！")
                        else:
                            st.warning(f"在相似度 > {sim_threshold}% 的条件下，未找到已上市的相似药物。你可以尝试调低相似度阈值。")
            else:
                st.warning("无法从 PDB 数据库获取该配体的 SMILES 序列，分析中止。")

except Exception as e:
    st.error(f"运行出错: {e}")
