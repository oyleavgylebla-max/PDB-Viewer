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
        ions = [
            'MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 
            'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL',
            'BR', 'I', 'F', 'IOD', 'FLC', 'NO3', 'NH4', 'ACT', 'FMT', 
            'EDO', 'GOL', 'PEG', 'DTT', 'BME'
        ]
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    return df

@st.cache_data
def get_smiles_from_pdb(ligand_id):
    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
        "Accept": "application/json"
    }
    try:
        pdbe_url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{ligand_id}"
        r = requests.get(pdbe_url, headers=headers, timeout=5)
        if r.status_code == 200:
            data = r.json()
            if ligand_id in data:
                smiles_list = data[ligand_id][0].get("smiles", [])
                for s in smiles_list:
                    if s.get("name") == "canonical": return s.get("value")
                if smiles_list: return smiles_list[0].get("value")
    except: pass 
    try:
        rcsb_url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
        r = requests.get(rcsb_url, headers=headers, timeout=5)
        if r.status_code == 200:
            data = r.json()
            for desc in data.get("rcsb_chem_comp_descriptor", []):
                if "SMILES" in desc.get("type", "").upper(): return desc.get("descriptor")
    except: pass
    return None

def search_chembl_drugs(smiles, similarity_threshold):
    """【Debug升级版】返回数据和具体的错误信息"""
    if not smiles or str(smiles).strip() == "":
        return [], "SMILES序列为空"
        
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json"
    
    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
        "Accept": "application/json"
    }
    
    try:
        # 将超时时间放宽到了 25 秒
        r = requests.get(url, headers=headers, timeout=25)
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            drugs = []
            for m in mols:
                # 放宽 max_phase 的判断条件，防止数据格式异常
                phase = m.get('max_phase', 0)
                if phase and float(phase) >= 4.0 and m.get('pref_name'):
                    drugs.append({
                        "药物名称 (Drug)": m.get('pref_name'),
                        "ChEMBL ID": m.get('molecule_chembl_id'),
                        "相似度 (%)": m.get('similarity', 'N/A'),
                        "分子量": m.get('molecule_properties', {}).get('full_mwt', 'N/A')
                    })
            return drugs, "Success"
        else:
            return [], f"API 拒绝访问，服务器返回状态码: {r.status_code}"
    except requests.exceptions.Timeout:
        return [], "请求超时！ChEMBL 数据库计算太慢，超过了25秒限制。"
    except Exception as e:
        return [], f"发生未知代码异常: {str(e)}"

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
                <div style="border: 1px solid #eee; border-radius: 10px; padding: 12px; margin-bottom: 20px; text-align: center; background: white; min-height: 280px;">
                    <strong style="color: #007bff; font-size: 1.1em;">{pdb_id}</strong>
                    {img_html}
                    <div style="font-size: 0.85em; color: #333; font-weight: bold; margin-bottom: 5px;">{target_id if target_id != 'ZZZ' else '-'}</div>
                    <div style="font-size: 0.75em; color: #666; height: 45px; overflow: hidden; line-height: 1.2;">{desc}</div>
                </div>
                """
                st.markdown(card_html, unsafe_allow_html=True)

    # ==========================================
    # 模式 B：单体详细查看 & AIDD 靶点发现
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
                st.image(pdbe_url, caption=f"配体 2D 结构图: {target_id}", width=250)
            with st.expander("查看完整描述"):
                st.write(info['Description (描述)'])
                st.write(f"文献出处: {info['Publication (文章出处)']}")

        with col2:
            st.subheader("3D 空间构象")
            view = py3Dmol.view(query=f'pdb:{selected_pdb.lower()}', width=800, height=500)
            view.setStyle({'cartoon': {'color': 'spectrum'}}) 
            view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
            view.zoomTo()
            showmol(view, height=500, width=800)

        # -----------------------------------------
        # 🚀 AIDD 分析引擎 (开启 Debug 显示)
        # -----------------------------------------
        if target_id != "ZZZ":
            st.divider()
            st.subheader("💊 AIDD 靶点与药物重定位分析")
            st.write("系统正在尝试寻找与当前 RNA 天然配体相似的 **FDA 已上市药物**...")
            
            auto_smiles = get_smiles_from_pdb(target_id)
            if not auto_smiles:
                st.warning(f"⚠️ 云端 API 被拦截或 {target_id} 暂无官方 SMILES。")
            
            smiles = st.text_input("🧬 配体 SMILES 序列 (自动抓取/支持手动覆盖):", value=auto_smiles if auto_smiles else "")
            
            if smiles:
                c1, c2 = st.columns([2, 1])
                with c1:
                    sim_threshold = st.slider("Tanimoto 结构相似度阈值 (%)", 70, 100, 70, 5)
                with c2:
                    st.write("")
                    st.write("")
                    search_btn = st.button("🚀 开始跨靶点搜索", use_container_width=True)
                
                if search_btn:
                    if not smiles.strip():
                        st.error("❌ SMILES 序列不能为空，请输入有效代码后再搜索！")
                    else:
                        with st.spinner("正在突破封锁，深入检索 ChEMBL 数据库 (最多可能需要 25 秒)..."):
                            # 获取数据和调试信息
                            fda_drugs, msg = search_chembl_drugs(smiles, sim_threshold)
                            
                            if fda_drugs:
                                st.success(f"🎉 找到 {len(fda_drugs)} 个相似的 FDA 上市药物！")
                                st.dataframe(pd.DataFrame(fda_drugs), use_container_width=True)
                                st.info("💡 **结论建议:** 如果这些老药原本属于其他疗法，它们可能通过相似的结合骨架与该 RNA 靶点结合。")
                            elif msg == "Success":
                                # 网络通了，但是真的没搜到这个药
                                st.warning("🔍 检索成功，但在当前阈值下，没有匹配到任何上市的 FDA 药物。")
                            else:
                                # 彻底失败，打印鲜红的错误信息
                                st.error(f"🚨 获取数据失败！详细错误原因：{msg}")
            else:
                st.info("👆 请在上方输入框手动填入该配体的 SMILES 序列后，开始匹配上市药物。")

except Exception as e:
    st.error(f"网页运行中发生错误: {e}")
