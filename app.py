import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol

# 1. 网页配置
st.set_page_config(page_title="RNA 结构精细化分类平台", layout="wide")
st.title("🧬 RNA 结构靶点分类与可视化系统")

@st.cache_data
def load_and_process_data():
    df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    
    # --- 核心升级：精细化分类引擎 ---
    def categorize(desc):
        desc_lower = str(desc).lower()
        # 1. 核糖开关
        if 'riboswitch' in desc_lower:
            return "核糖开关 (Riboswitch)"
        # 2. 适配体
        elif 'aptamer' in desc_lower:
            return "适配体 (Aptamer)"
        # 3. G-四联体 (新增)
        elif any(word in desc_lower for word in ['quadruplex', 'g-4', 'g4']):
            return "G-四联体 (G-quadruplex)"
        # 4. 核糖体 (从核酶中独立出来)
        elif any(word in desc_lower for word in ['ribosomal', 'ribosome', 'rrna']):
            return "核糖体 (rRNA)"
        # 5. 核酶
        elif 'ribozyme' in desc_lower:
            return "核酶 (Ribozyme)"
        # 6. 特殊结构基元
        elif any(word in desc_lower for word in ['ires', 'hairpin', 'stem-loop', 'pseudoknot']):
            return "特殊结构基元 (Special/Motifs)"
        else:
            return "其他 RNA (Others)"
            
    df['Category'] = df['Description (描述)'].apply(categorize)
    
    # 为每一行预先提取主配体 ID，用于后续的“相似性聚类排序”
    def get_sort_id(ligand_text):
        if ligand_text == "No ligands" or str(ligand_text) == "nan": return "ZZZ" # 没配体的排最后
        all_l = str(ligand_text).split(' | ')
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'BR', 'I', 'F', 'EDO', 'GOL']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    return df

try:
    df = load_and_process_data()
    
    # 2. 侧边栏菜单
    st.sidebar.header("⚙️ 显示模式")
    view_mode = st.sidebar.radio("模式切换", ["📊 同类全局对比 (画廊)", "🔍 单体详细查看 (3D)"])
    
    st.sidebar.divider()
    st.sidebar.header("🔍 分类筛选")
    category_list = ["全部 (All)"] + sorted(df['Category'].unique().tolist())
    selected_cat = st.sidebar.selectbox("选择 RNA 类型", category_list)
    
    # 筛选并按照配体 ID 排序（实现相似分子聚在一块）
    if selected_cat == "全部 (All)":
        filtered_df = df.sort_values(by=['MainLigandID', 'PDB ID'])
    else:
        filtered_df = df[df['Category'] == selected_cat].sort_values(by=['MainLigandID', 'PDB ID'])

    # ==========================================
    # 模式 A：同类全局对比 (画廊模式)
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
                <div style="border: 1px solid #eee; border-radius: 10px; padding: 12px; margin-bottom: 20px; text-align: center; background: white; box-shadow: 2px 2px 5px rgba(0,0,0,0.05);">
                    <strong style="color: #007bff; font-size: 1.1em;">{pdb_id}</strong>
                    {img_html}
                    <div style="font-size: 0.85em; color: #333; font-weight: bold; margin-bottom: 5px;">{target_id if target_id != 'ZZZ' else '-'}</div>
                    <div style="font-size: 0.75em; color: #666; height: 45px; overflow: hidden; line-height: 1.2;">{desc}</div>
                </div>
                """
                st.markdown(card_html, unsafe_allow_html=True)

    # ==========================================
    # 模式 B：单体详细查看 (3D 模式)
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
                st.image(pdbe_url, caption=f"配体 2D 结构: {target_id}", width=250)
            
            with st.expander("查看完整描述"):
                st.write(info['Description (描述)'])
                st.write(f"文献: {info['Publication (文章出处)']}")

        with col2:
            st.subheader("3D 空间构象")
            view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=550)
            view.setStyle({'cartoon': {'color': 'spectrum'}}) 
            view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
            view.zoomTo()
            showmol(view, height=550, width=800)

except Exception as e:
    st.error(f"运行出错: {e}")
