import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol

# 1. 网页配置
st.set_page_config(page_title="PDB RNA 结构可视化平台", layout="wide")
st.title("🧬 PDB 结构分类与可视化分析")

@st.cache_data
def load_and_process_data():
    df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    
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

# 辅助函数：提取真正的小分子ID并生成图片链接
def get_ligand_img_url(ligand_text):
    if ligand_text == "No ligands" or str(ligand_text) == "nan":
        return None
    
    all_ligands = str(ligand_text).split(' | ')
    
    # 🚀 扩充版黑名单：涵盖金属离子、卤素、常见结晶辅料和缓冲液成分
    ions = [
        'MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 
        'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', # 金属与碱土
        'BR', 'I', 'F', 'IOD', 'FLC', # 卤素离子
        'NO3', 'NH4', 'ACT', 'FMT', 'EDO', 'GOL', 'PEG', 'DTT', 'BME', 'GOL' # 结晶溶剂/缓冲液
    ]
    
    target_id = None
    for l in all_ligands:
        lid = l.strip().split(' ')[0].upper()
        if lid not in ions:
            target_id = lid
            break
            
    # 如果找了一圈全是黑名单里的东西，说明可能真的是个只结合金属的结构
    if not target_id:
        target_id = all_ligands[0].strip().split(' ')[0].upper()
        
    return target_id

try:
    df = load_and_process_data()
    
    # 2. 侧边栏菜单设计
    st.sidebar.header("⚙️ 查看设置")
    
    view_mode = st.sidebar.radio("选择显示模式", ["🔍 单体详细查看 (3D)", "📊 同类全局对比 (2D画廊)"])
    
    st.sidebar.divider()
    st.sidebar.header("🔍 数据筛选")
    
    category_list = ["全部 (All)"] + sorted(df['Category'].unique().tolist())
    selected_cat = st.sidebar.selectbox("选择结构分类", category_list)
    
    if selected_cat == "全部 (All)":
        filtered_df = df
    else:
        filtered_df = df[df['Category'] == selected_cat]

    # ==========================================
    # 模式 A：同类全局对比 (画廊模式)
    # ==========================================
    if view_mode == "📊 同类全局对比 (2D画廊)":
        st.subheader(f"当前分类: {selected_cat} (共 {len(filtered_df)} 个结构)")
        
        with st.expander("📝 展开查看完整数据表格"):
            st.dataframe(filtered_df[['PDB ID', 'Description (描述)', 'Ligands (对应小分子)']], use_container_width=True)
        
        st.write("### 🖼️ 核心配体对比画廊")
        
        cols_per_row = 4
        columns = st.columns(cols_per_row)
        
        for idx, row in filtered_df.reset_index().iterrows():
            col_idx = idx % cols_per_row
            
            pdb_id = row['PDB ID']
            desc = row['Description (描述)']
            ligand_text = row['Ligands (对应小分子)']
            
            target_id = get_ligand_img_url(ligand_text)
            
            with columns[col_idx]:
                if target_id:
                    img_url = f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{target_id}_400.svg"
                    img_html = f'<img src="{img_url}" style="width:100%; height:150px; object-fit:contain; margin-bottom:10px;">'
                else:
                    img_html = '<div style="height:150px; display:flex; align-items:center; justify-content:center; color:gray; background:#f0f2f6; border-radius:5px; margin-bottom:10px;">无配体</div>'
                
                card_html = f"""
                <div style="border: 1px solid #e6e6e6; border-radius: 8px; padding: 10px; margin-bottom: 20px; text-align: center; background-color: white;">
                    <h4 style="margin-top:0; color: #1f77b4;">{pdb_id}</h4>
                    {img_html}
                    <div style="font-size: 0.8em; color: #555; height: 60px; overflow: hidden; text-overflow: ellipsis;">
                        <b>{target_id if target_id else 'N/A'}</b><br/>
                        {desc}
                    </div>
                </div>
                """
                st.markdown(card_html, unsafe_allow_html=True)


    # ==========================================
    # 模式 B：单体详细查看 (3D 模式)
    # ==========================================
    else:
        selected_pdb = st.sidebar.selectbox("选择 PDB ID", filtered_df['PDB ID'].tolist())
        info = df[df['PDB ID'] == selected_pdb].iloc[0]
        
        st.sidebar.divider()
        st.sidebar.subheader("🧪 小分子 2D 结构式")
        
        ligand_text = str(info['Ligands (对应小分子)'])
        target_id = get_ligand_img_url(ligand_text)
        
        if target_id:
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

        col1, col2 = st.columns([1, 2])

        with col1:
            st.subheader("详细信息")
            st.success(f"**结构分类:** {info['Category']}")
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
