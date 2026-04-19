import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol

# 1. 网页配置
st.set_page_config(page_title="PDB Riboswitch 可视化平台", layout="wide")
st.title("🧬 PDB 结构分类与可视化分析")

@st.cache_data
def load_and_process_data():
    # 读取表格
    df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    
    # --- 新增分类逻辑 ---
    def categorize(desc):
    desc_lower = str(desc).lower()
    if 'riboswitch' in desc_lower:
        return "核糖开关 (Riboswitch)"
    elif 'aptamer' in desc_lower:
        return "适配体 (Aptamer)"
    elif any(word in desc_lower for word in ['ribozyme', 'ribosomal', 'ribosome']):
        return "催化核酶 (Ribozyme/Ribosomal)"
    elif any(word in desc_lower for word in ['ires', 'hairpin', 'stem-loop']):
        return "特殊调控/结构基元 (Special/Motifs)"
    else:
        return "其他 RNA 结构 (Other RNA)"
            
    df['Category'] = df['Description (描述)'].apply(categorize)
    return df

try:
    df = load_and_process_data()
    
    # 2. 侧边栏分类筛选
    st.sidebar.header("🔍 数据筛选")
    
    # 先选大类
    category_list = ["全部 (All)"] + sorted(df['Category'].unique().tolist())
    selected_cat = st.sidebar.selectbox("选择结构分类", category_list)
    
    # 根据大类筛选 ID 列表
    if selected_cat == "全部 (All)":
        filtered_df = df
    else:
        filtered_df = df[df['Category'] == selected_cat]
        
    selected_pdb = st.sidebar.selectbox("选择 PDB ID", filtered_df['PDB ID'].tolist())
    
    # 获取选中行数据
    info = df[df['PDB ID'] == selected_pdb].iloc[0]
    
    # 3. 侧边栏 2D 结构图 (沿用之前的 PDBe 稳定版接口)
    st.sidebar.divider()
    st.sidebar.subheader("🧪 小分子 2D 结构式")
    
    ligand_text = str(info['Ligands (对应小分子)'])
    if ligand_text != "No ligands" and ligand_text != "nan":
        all_ligands = ligand_text.split(' | ')
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA']
        
        target_id = None
        for l in all_ligands:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions:
                target_id = lid
                break
        if not target_id:
            target_id = all_ligands[0].strip().split(' ')[0].upper()

        st.sidebar.markdown(f"**主要配体 ID:** `{target_id}`")
        pdbe_img_url = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg"
        html_code = f"""
        <div style="background-color: white; border: 1px solid #ddd; border-radius: 8px; padding: 10px; text-align: center;">
            <img src="{pdbe_img_url}" alt="Ligand {target_id}" style="max-width: 100%; max-height: 300px; object-fit: contain;">
        </div>
        """
        st.sidebar.markdown(html_code, unsafe_allow_html=True)
    else:
        st.sidebar.info("该结构不含小分子配体")

    # 4. 主界面展示
    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("详细信息")
        st.success(f"**分类:** {info['Category']}") # 突出显示分类
        st.info(f"**PDB ID:** {selected_pdb}")
        st.write(f"**描述:** {info['Description (描述)']}")
        st.write(f"**文献:** {info['Publication (文章出处)']}")
        st.write(f"**配体/离子:** {info['Ligands (对应小分子)']}")

    with col2:
        st.subheader("3D 交互视图")
        view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=500)
        view.setStyle({'cartoon': {'color': 'spectrum'}}) 
        view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
        view.zoomTo()
        showmol(view, height=500, width=800)

except Exception as e:
    st.error(f"发生错误: {e}")
