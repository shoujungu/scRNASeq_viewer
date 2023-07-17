import PySimpleGUI as sg
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import anndata as ad
import numpy as np
import os

sg.theme('SystemDefaultForReal')
l_reset=[]

#-----------------------------------------------------
#layout
layout=[
[sg.T('')], 
[sg.T('Input File:', font=8, key='-fname-'), 
    sg.Input(key='-fin-', enable_events=True, visible=False),
    sg.FileBrowse()],

[sg.T('        Output Folder:', font=8), 
    sg.FolderBrowse(key='-fdout-'),
    sg.Button('Reset', size=(6,1))],

[sg.T('Gene:', font=8), 
    sg.Input(size=(15, 1), key='-gene-'), 
    sg.Button('Submit')],

[sg.Button('Exit', size=(6, 1))],
[sg.HSeparator(color='#888888')],

[sg.Listbox(values=[], enable_events=True, size=(20, 700), key='-cell-', font=6), 
    sg.Image(key='-anno-', size=(860, 700)), 
    sg.VSeparator(color='#000000'), 
    sg.Image(key='-exp-', size=(950, 700))]
]

#window
window=sg.Window('UMAP Viewer', layout, size=(1900,1000), element_justification="center")

#-----------------------------------------------------
def plt_clus(ada, cell, f_out, check, col='cell', sz=(2.5, 2.5), cmap=['#ff0000', '#999999']):
    s=round(4000/ada.shape[0], 1)
    s=s if s>0.5 else 0.5
    s=s if s<2 else 2
    #prep
    adai=ada.copy()
    adai.obs['clus']=adai.obs[col].apply(lambda x: cell if x==cell else 'Others') 
    df=pd.DataFrame(adai.obsm['X_umap'], columns=['x', 'y'], index=adai.obs.index)
    df['clus']=pd.Categorical(adai.obs['clus'], categories=[cell, 'Others'], ordered=True)
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue='clus', palette=cmap, s=s, linewidth=0, alpha=0.8)
    #adjust
    plt.axis('off')
    ax.set_title(cell, weight='semibold')
    ax.get_legend().remove()
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    check=True
    return check


def plt_exp(ada, gene, f_out, check, sz=(2.8, 2.5), cmap='BuPu', vmin=-0.5):
    s=round(4000/ada.shape[0], 1)
    s=s if s>0.5 else 0.5
    s=s if s<2 else 2
    #prep
    df_xy=pd.DataFrame(ada.obsm['X_umap'], columns=['x', 'y'], index=ada.obs.index)
    df=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var).to_df()
    df=df.merge(df_xy, left_index=True, right_index=True)
    #norm cmap
    vmax=df[gene].quantile(1)
    norm=plt.Normalize(vmin, vmax)
    l_tick=np.arange(-0.5, vmax, 0.5).tolist()
    #plot  (left, bottom, width, height)
    fig=plt.figure(figsize=sz)
    ax=fig.add_axes([0, 0.05, 0.85, 0.8])
    ax=sns.scatterplot(x='x', y='y', data=df, hue=gene, palette=cmap, s=s, linewidth=0, alpha=0.8, hue_norm=norm)
    #adjust
    plt.axis('off')
    ax.set_title(gene, weight='semibold')
    ax.get_legend().remove()
    #color bar (https://stackoverflow.com/questions/62884183/trying-to-add-a-colorbar-to-a-seaborn-scatterplot)
    ax_cbar=fig.add_axes([0.85, 0.3, 0.02, 0.4])
    sm=plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb=fig.colorbar(sm, cax=ax_cbar, ticks=l_tick)
    cb.ax.tick_params(labelsize=5)
    #save
    plt.savefig(f_out, dpi=300)
    plt.close()
    check=True
    return check


def e_fin(fin): #update listbox
    try:
        name=Path(fin).stem
        ada=sc.read(fin)
        l_cell=ada.obs['cell'].astype('str').unique().tolist()
        l_cell.sort()
        window['-cell-'].update(l_cell)
        window['-fname-'].update(f'Input File: {name}')
    except Exception as e:
        window['-cell-'].update([])
        print(e)
    return


def e_reset(fd_out):  #delete all png
    l_png=list(Path(fd_out).glob('*.png'))
    for f in l_png: os.remove(f)
    return


def e_cell(fin, fd_out, cell): #plot cell cluster
    check=False
    try:
        #prep
        ada=sc.read(fin)
        cell=cell[0]
        #plot
        f_out=f'{fd_out}/{cell}.png'
        check=plt_clus(ada, cell, f_out, check)
    except Exception as e:  print(e)
    #update image
    if check: window['-anno-'].update(filename=f_out)
    return


def e_submit(fin, fd_out, gene):  #plot gene exp
    check=False
    try:
        ada=sc.read(fin)
        #plot
        f_out=f'{fd_out}/{gene}.png'
        check=plt_exp(ada, gene, f_out, check)
    except Exception as e:  print(e)
    #update image
    if check: window['-exp-'].update(filename=f_out)
    return

#-----------------------------------------------------
while True:
    event, d_val=window.read()
    if event==sg.WIN_CLOSED or event=="Exit": break

    #update listbox when select input
    elif event=='-fin-': e_fin(d_val['-fin-'])
    
    #delete all png files in the output folder
    elif event=='Reset': e_reset(d_val['-fdout-'])

    #plot cell cluster
    elif event=='-cell-': e_cell(d_val['-fin-'], d_val['-fdout-'], d_val['-cell-'])
    
    #plot exp
    elif event=='Submit': e_submit(d_val['-fin-'], d_val['-fdout-'], d_val['-gene-'])
