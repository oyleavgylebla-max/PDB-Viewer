import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol

# 1. 设置网页标题
st.set_page_config(page_title="PDB 结构与小分子可视化", layout="wide")
st.title("🧬 PDB 结构可视化与化学结构分析")

@st.cache_data
def load_data():
    # 改成这样：
    return pd.read_excel("PDB_Dataset_Info_Full.xlsx")

try:
    df = load_data()
    
    st.sidebar.header("🔍 结构筛选")
    selected_pdb = st.sidebar.selectbox("选择 PDB ID", df['PDB ID'].tolist())
    
    info = df[df['PDB ID'] == selected_pdb].iloc[0]
    
    st.sidebar.divider()
    st.sidebar.subheader("🧪 小分子 2D 结构式")
    
    ligand_text = str(info['Ligands (对应小分子)'])
    
    if ligand_text != "No ligands" and ligand_text != "nan":
        all_ligands = ligand_text.split(' | ')
        # 黑名单：过滤常见离子
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

        # -------------------------------------------------------------
        # 🚀 绝杀技：切换到欧洲 PDB (PDBe) 的开源图片服务器！
        # 他们允许任何网页直接调用高清 SVG 矢量图，再也不会有防盗链拦截！
        # -------------------------------------------------------------
        pdbe_img_url = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg"
        
        html_code = f"""
        <div style="background-color: white; border: 1px solid #ddd; border-radius: 8px; padding: 10px; text-align: center;">
            <img src="{pdbe_img_url}" 
                 alt="Ligand {target_id}"
                 style="max-width: 100%; max-height: 300px; object-fit: contain;">
        </div>
        """
        st.sidebar.markdown(html_code, unsafe_allow_html=True)
            
    else:
        st.sidebar.info("该结构不含任何小分子配体。")

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("详细信息")
        st.info(f"**PDB ID:** {selected_pdb}")
        st.write(f"**描述:** {info['Description (描述)']}")
        st.write(f"**文献:** {info['Publication (文章出处)']}")
        st.write(f"**包含配体/离子:** {info['Ligands (对应小分子)']}")

    with col2:
        st.subheader("3D 交互视图 (可鼠标拖拽、滚轮缩放)")
        view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=500)
        view.setStyle({'cartoon': {'color': 'spectrum'}}) 
        view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
        view.zoomTo()
        showmol(view, height=500, width=800)

except Exception as e:
    st.error(f"发生错误: {e}")