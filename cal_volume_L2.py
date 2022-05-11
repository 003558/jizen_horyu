
# coding: utf-8
from winreg import QueryInfoKey
import yaml
import glob
import shutil
from datetime import datetime, timedelta
import numpy as np
import subprocess        # import文なので次以降の例では省略します
import pandas as pd
import csv
import datetime
from scipy import signal, interpolate
import matplotlib.pyplot as plt
import os
from io import StringIO
import math
from ctypes import *
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d
# グラフ用のフォントの設定
from matplotlib.font_manager import FontProperties
# x軸の日付表示
import matplotlib.dates as mdates
# グラフのサイズ変更
from matplotlib import gridspec
# グラフ用のフォントの設定
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Hiragino Maru Gothic Pro', 'Yu Gothic', 'Meirio', 'Takao', 'IPAexGothic', 'IPAPGothic', 'VL PGothic', 'Noto Sans CJK JP']
fp1 = FontProperties(fname="C:\Windows\Fonts\msgothic.ttc", size=15)
fp2 = FontProperties(fname="C:\Windows\Fonts\msgothic.ttc", size=12)
# 実行時の警告を表示しない
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import warnings
warnings.simplefilter('ignore')

def read_ame(file, extra):
    area_9 = 33.9
    area_19 = 87.1
    area_20 = 164.5
    area_41 = 66.2
    f = open(file, 'r', encoding='UTF-8')
    datalist = f.readlines()
    length = int(datalist[0].split()[0])
    if length % 8 == 0:
        row = int(length / 8)
    else:
        row = int(length / 8)
    i = 1
    df_r = pd.DataFrame()
    print(length)
    for data in range(length):
        ryuiki = int(datalist[i].split()[0])
        i += 1
        if ryuiki == 9 or ryuiki == 19 or ryuiki == 20 or ryuiki == 41:
            r_list = [0.0 for x in range(extra)]
            for n in range(row):
                r_list += [float(x) for x in datalist[i].split()]
                i += 1
            df_r[str(ryuiki)] = r_list
        else:
            i += row
        if data == 45:
            break
    f.close()
    df_r["r"] = (area_9*df_r["9"]+area_19*df_r["19"]+area_20*df_r["20"]+area_41*df_r["41"])/(area_9+area_19+area_20+area_41)

    return df_r

def time_naiso(df, term, term_moto, ini_time):
    # 
    dt = timedelta(minutes=term)
    num_new = int(len(df)*term_moto/term-1)
    x = df['Qin'].values
    y = df['flg_2ji'].values
    t = pd.to_datetime(df['date'].values)
    t_new = [ini_time + x * dt for x in range(num_new)]

    # interp関数はdatetime型を受け付けない
    # 時刻をdatetime型からunix時間（float）に変換する
    t_unix = [x.timestamp() for x in t]
    t_new_unix = [x.timestamp() for x in t_new]
    #t_new = [datetime.datetime.fromtimestamp(x) for x in t_new_unix]

    # 方法1: numpy で補間
    x_numpy = np.interp(t_new_unix, t_unix, x)
    y_numpy = np.interp(t_new_unix, t_unix, y)

    df_new = pd.DataFrame()
    df_new['date'] = t_new
    df_new['Qin']  = x_numpy
    df_new['flg_2ji']  = y_numpy
    
    return df_new

class dam():
    def __init__(self):
        self.outf      = 'Output' # 出力フォルダ
        
        #各種データのファイル名設定
        self.path_gensoku = './amagase_dam/gensoku.csv'    #放流の原則テーブル!!!
        self.path_HV = './amagase_dam/HV.csv'      #HVテーブル!!!
        self.path_HQ = './amagase_dam/HQ.csv'      #HQテーブル!!!
        self.path_Qin = './amagase_dam/Qin.csv'

        self.path_res= './amagase_dam/res/res.csv'                #利水計算結果
        self.path_grhp='./amagase_dam/res/praph.png'              #利水計算結果グラフ
        
        self.gensoku    = pd.read_csv(self.path_gensoku)
        self.HV_table   = pd.read_csv(self.path_HV)
        self.HQ_table   = pd.read_csv(self.path_HQ)

        #HV補間関数作成
        HVH_list          = self.HV_table.iloc[0:,0].values.tolist()
        HVV_list          = (self.HV_table.iloc[0:,1]).values.tolist()
        self.interpolation_HV  = interp1d(HVH_list, HVV_list, kind="linear")
        self.interpolation_VH  = interp1d(HVV_list, HVH_list, kind="linear")
        
        #HQ補間関数作成
        HQH_list          = self.HQ_table.iloc[0:,0].values.tolist()
        HQQ_list          = (self.HQ_table.iloc[0:,1]).values.tolist()
        self.interpolation_HQ  = interp1d(HQH_list, HQQ_list, kind="linear")
        self.interpolation_QH  = interp1d(HQQ_list, HQH_list, kind="linear")
        
        #tadasigaki
        self.Ht = 76.0
        self.HtV = self.interpolation_HV(self.Ht)
        
        #サーチャージ
        self.Hser = 78.5
        self.HserV = self.interpolation_HV(self.Hser)
        
        #制限水位
        Hs = 72.0
        self.HsV = self.interpolation_HV(Hs)
        
    """■笂ｵ諸々の変数を初期化して洪水イベントを取得"""
    def _reset(self, Hini, extra):                                                                                    # 洪水イベントを抽出
        self.steps     = 0
        self.Qin_df_tmp= pd.read_csv(self.path_Qin)
        ini_time = pd.to_datetime(self.Qin_df_tmp['date'].values)[0] - datetime.timedelta(hours=extra)
        self.Qin_peak  = max(self.Qin_df_tmp['Qin'].values)

        #add extra
        extra_list = []
        for i in range(extra):
            extra_list.append([self.Qin_df_tmp['date'].values[0], self.Qin_df_tmp['Qin'].iloc[0], self.Qin_df_tmp['flg_2ji'].iloc[0]])
        extra_df = pd.DataFrame(extra_list, columns=['date','Qin','flg_2ji'])
        self.Qin_df_tmp = pd.concat([extra_df,self.Qin_df_tmp])

        #Qin naiso
        self.Qin_df   = time_naiso(self.Qin_df_tmp, 15, 60, ini_time)
        self.max_steps = len(self.Qin_df)
        self.Qin_df.to_csv("./aaa.csv")
        self.steps    = 0
        self.Qin      = self.Qin_df['Qin'].values[0]
        self.V        = self.interpolation_HV(Hini)
        self.Qout     = 0.1
        self.Qout_pre = 0.1
        self.Qin_pre  = self.Qin
        self.H        = Hini
        self.H_pre    = self.H
        self.V_pre    = self.V
        self.Vini     = self.V
        self.peak_flg = 0
        self.flg_2ji_pre  = 0
        self.flg_2ji  = 0
        self.kouki    = 0 
        self.y_step    = -999  #予備放流開始前後
        self.j_step    = -999  #事前放流開始前後
        self.r_step    = -999  #基準雨量超過前後
        self.junbi_step= 0    
        self.tujo_flg  = 0 

        #output用
        self.output_Qin_o  = []
        self.output_Qout_o = []
        self.output_H_o    = []
        self.output_Qout_p = []
        self.output_H_p    = []
        self.output_time   = []
        self.output_rain   = []
        self.output_sum_r_y= []
        self.output_sum_r_c= []


    def _step(self, step, df_r):# 初回1ステップ実行
        
        #現時刻のダム諸量
        self.H_pre = self.H  
        self.V_pre =  self.interpolation_HV(self.H)
        self.Qout_pre = self.Qout
        self.Qin_pre  = self.Qin
        self.flg_2ji_pre  = self.flg_2ji
        self.flg_2ji  = self.Qin_df['flg_2ji'].values[step]

        #ruika-uryo
        self.sum_rain_y = 0.0
        self.sum_rain_c = 0.0
        for i in range(1, 40):
            if int(step/4)+i < len(df_r):
                self.sum_rain_y += (df_r['r'].iloc[int(step/4)+i])
            else:
                self.sum_rain_y += 0.0
            if int(step/4)-i > 0 and int(step/4)-i < len(df_r):
                self.sum_rain_c += (df_r['r'].iloc[int(step/4)-i])
            else:
                self.sum_rain_c += 0.0

        #   Qout設定
        self.Qin = self.Qin_df['Qin'].values[step]
        if self.tujo_flg == 1:
            if self.Qout < 1140.0:
                Qout_gensoku_max = self.gensoku[(self.gensoku['Qout_d']<=self.Qout) & (self.gensoku['Qout_u']>self.Qout)]['dQout'].values[0]
            else:
                Qout_gensoku_max = 0.0
            if int(self.flg_2ji) == 1:
                self.Qout = 250.0
            elif int(self.flg_2ji) == 0 and int(self.flg_2ji_pre) == 1:
                self.Qout = min(self.Qin, 1140.0)
                self.kouki = 1
            else:
                if self.kouki == 0:
                    self.Qout = min(self.Qin, Qout_gensoku_max+self.Qout, self.interpolation_HQ(self.H), 1140.0)
                else:
                    self.Qout = min(self.Qin, 1140.0)
            if self.H >= self.Ht:
                self.Qout = self.Qin
                self.H = self.Ht
            #次ステップ水位算定
            self.V     = self.V_pre + (self.Qin - self.Qout) * 900
            self.H     =  self.interpolation_VH(self.V)
        else:
            if self.y_step<0 and self.r_step<0 and self.j_step<0:
                self.Qout = self.Qin
                if self.sum_rain_y > 150.0:
                    self.r_step = self.steps
                    self.r_step_x = [pd.to_datetime(self.Qin_df['date'].values[step]).strftime('%m/%d %H:%M'), pd.to_datetime(self.Qin_df['date'].values[step]).strftime('%m/%d %H:%M')]
                    self.r_step_y = [0, 10000]
            elif self.y_step<0 and not(self.r_step<0) and self.j_step<0:
                self.Qout = self.Qin
                self.junbi_step += 1
                if self.junbi_step == 6*4+1:
                    self.y_step = self.steps
                    self.y_step_x = [pd.to_datetime(self.Qin_df['date'].values[step]).strftime('%m/%d %H:%M'), pd.to_datetime(self.Qin_df['date'].values[step]).strftime('%m/%d %H:%M')]
                    self.y_step_y = [0, 10000]
            elif not(self.y_step<0) and not(self.r_step<0) and self.j_step<0:
                if self.Qout < 1140.0:
                    Qout_gensoku_max = self.gensoku[(self.gensoku['Qout_d']<=self.Qout) & (self.gensoku['Qout_u']>self.Qout)]['dQout'].values[0]
                    self.Qout = min(Qout_gensoku_max+self.Qout, self.interpolation_HQ(self.H), 1140.0)
                else:
                    self.Qout = 1140.0
                #次ステップ水位算定
                self.V     = self.V_pre + (self.Qin - self.Qout) * 900
                self.H     =  self.interpolation_VH(self.V)
                if self.H < 58.0:
                    self.j_step = self.steps
                    self.j_step_x = [pd.to_datetime(self.Qin_df['date'].values[step]).strftime('%m/%d %H:%M'), pd.to_datetime(self.Qin_df['date'].values[step]).strftime('%m/%d %H:%M')]
                    self.j_step_y = [0, 10000]
                if self.Qout_pre > self.Qout:
                    self.tujo_flg = 1
            elif not(self.y_step<0) and not(self.r_step<0) and not(self.j_step<0):
                if self.Qout < 1140.0:
                    Qout_gensoku_max = self.gensoku[(self.gensoku['Qout_d']<=self.Qout) & (self.gensoku['Qout_u']>self.Qout)]['dQout'].values[0]
                    self.Qout = min(Qout_gensoku_max+self.Qout, self.interpolation_HQ(self.H), 1140.0)
                else:
                    self.Qout = 1140.0
                #次ステップ水位算定
                self.V     = max(self.V_pre + (self.Qin - self.Qout) * 900, self.HV_table.iloc[0,1])
                self.H     =  self.interpolation_VH(self.V)
                if self.Qin > self.Qout:
                    self.tujo_flg = 1
                if self.H < 45.1:
                    self.H = 45.1
                    self.tujo_flg = 1
            
                
        #計算諸量
        if self.Qin_peak == self.Qin:
            self.peak_flg = 1
        else:
            self.peak_flg = 1
        
        #output用
        self.output_Qin_o.append(self.Qin)
        self.output_Qout_p.append(self.Qout)
        self.output_H_p.append(self.H)
        self.output_sum_r_y.append(self.sum_rain_y)
        self.output_sum_r_c.append(self.sum_rain_c)
        self.output_time.append(pd.to_datetime(self.Qin_df['date'].values[step]).strftime('%m/%d %H:%M'))
        if step <= len(df_r)*4-1:
            self.output_rain.append(df_r['r'].iloc[int(step/4)])
        else:
            self.output_rain.append(0.0)

        #if (self.steps>=self.max_steps-1) or ((self.Qin_pre <= self.Qout_pre) and (self.Hs <= self.H)) or ((self.Qin_pre <= self.Qout_pre) and (self.peak_flg == 1)):
        #    self.done = 1
        


    def _loop(self, df_r, name):
        for step in range(self.max_steps):
            self._step(step, df_r)
        print("------------")


    def _graph(self,title):
        self.obs_all       = []
        self.obs_all.append([self.output_time,
            self.output_H_p,
            self.output_Qin_o,
            self.output_Qout_p,
            self.output_rain,
            self.output_sum_r_y,
            self.output_sum_r_c])
        
        self.df_obs_all = pd.DataFrame(self.obs_all, columns=['date', 'H', 'QinDam', 'QoutDam', 'rain', 'sum_r_y', 'sum_r_c'])
        self.df_obs_all = self.df_obs_all.set_index('date')
        #out_csv
        df_csv = pd.DataFrame()
        df_csv['time'] = self.output_time
        df_csv['H'] = self.output_H_p
        df_csv['Qin'] = self.output_Qin_o
        df_csv['Qout'] = self.output_Qout_p
        df_csv['rain'] = self.output_rain
        df_csv.to_csv("check.csv")
                                
        #
        fig = plt.figure(figsize=(24.0, 16.0))
        grid = plt.GridSpec(4, 1, wspace=0.1, hspace=0.2)
        ax1 = plt.subplot(grid[2:4, 0])
        ln1=ax1.plot(self.df_obs_all.index.values[0], self.df_obs_all['QinDam'].values[0], label='流入量')
        ln1=ax1.plot(self.df_obs_all.index.values[0], self.df_obs_all['QoutDam'].values[0], label='放流量')
        ax2 = ax1.twinx()
        ln2=ax2.plot(self.df_obs_all.index.values[0], self.df_obs_all['H'].values[0],label='貯水位',color='black', zorder=1)
        ln2=ax2.plot([self.df_obs_all.index.values[0][0],self.df_obs_all.index.values[0][-1]], [self.Hser, self.Hser],color='red', linestyle="--")
        ln2=ax2.plot([self.df_obs_all.index.values[0][0],self.df_obs_all.index.values[0][-1]], [self.Ht, self.Ht],color='orange', linestyle="--")
        ln2=ax2.plot(self.r_step_x,self.r_step_y,color='green', linestyle="-")
        ln2=ax2.plot(self.y_step_x,self.y_step_y,color='green', linestyle="--")
        ln2=ax2.plot(self.j_step_x,self.j_step_y,color='green', linestyle="dashdot")
        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax1.set_ylim(0, 6000, 500)                                                            # y軸範囲の指定
        ax1.yaxis.set_ticks(np.arange(0, 6000+0.0001, 500))                                 # y軸範囲の指定
        ax1.xaxis.set_ticks(np.arange(0, len(self.output_Qout_p), step=48))
        ax2.set_ylim(44, 80, 3)                                                               # y軸範囲の指定
        ax2.yaxis.set_ticks(np.arange(44, 80+0.0001, 3))                                 # y軸範囲の指定
        #ax1.set_xlabel('time')
        ax1.set_ylabel(r'流量(㎥/s)', fontsize=18)
        ax1.grid(True)
        ax2.set_ylabel(r'水位(O.P.m)', fontsize=18)
        ax2.legend(h1+h2, l1+l2, bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=0, edgecolor='black', framealpha=1, facecolor='w', prop={"family":"MS Gothic"}, fontsize=18)
        ax1.axes.xaxis.set_visible(True)
        ax1.xaxis.grid(b=True,which='major',color = "gray", alpha = 0.8)
        ax1.xaxis.grid(b=True,which='minor',color = "lightgray", alpha = 0.8, linestyle = "--")
        ax1.yaxis.grid(b=True,which='major',color = "gray", alpha = 0.8)
        ax1.tick_params(labelsize=18)
        ax2.tick_params(labelsize=18)
        
        ax3 = plt.subplot(grid[0, 0])
        ln3 = ax3.bar(self.df_obs_all.index.values[0], self.df_obs_all['rain'].values[0], width=1.0, label='引き伸ばし雨量')
        ln3 = ax3.plot(self.r_step_x,self.r_step_y,color='green', linestyle="-")
        ln3 = ax3.plot(self.y_step_x,self.y_step_y,color='green', linestyle="--")
        ln3 = ax3.plot(self.j_step_x,self.j_step_y,color='green', linestyle="dashdot")
        ax3.set_ylim(0, 120, 20)                                                            # y軸範囲の指定
        ax3.yaxis.set_ticks(np.arange(0, 120+0.0001, 20))                                 # y軸範囲の指定
        ax3.axes.xaxis.set_visible(False)
        ax3.invert_yaxis()
        ax3.xaxis.set_ticks(np.arange(0, len(self.output_Qout_p), step=48))
        ax3.grid(True)
        ax3.set_ylabel(r'雨量(mm/hr)', fontsize=18)
        ax3.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=0, edgecolor='black', facecolor='w', framealpha=1, prop={"family":"MS Gothic"})
        ax3.tick_params(labelsize=18)

        ax4 = plt.subplot(grid[1, 0])
        ln4 = ax4.plot(self.df_obs_all.index.values[0], self.df_obs_all['sum_r_c'].values[0], label='実績累加39時間雨量')
        ln4 = ax4.plot(self.df_obs_all.index.values[0], self.df_obs_all['sum_r_y'].values[0], label='予測累加39時間雨量')
        ln4 = ax4.plot(self.r_step_x,self.r_step_y,color='green', linestyle="-")
        ln4 = ax4.plot(self.y_step_x,self.y_step_y,color='green', linestyle="--")
        ln4 = ax4.plot(self.j_step_x,self.j_step_y,color='green', linestyle="dashdot")
        ax4.set_ylim(0, 700, 100)                                                            # y軸範囲の指定
        ax4.yaxis.set_ticks(np.arange(0, 700+0.0001, 100))                                 # y軸範囲の指定
        ax4.axes.xaxis.set_visible(False)
        ax4.invert_yaxis()
        ax4.xaxis.set_ticks(np.arange(0, len(self.output_Qout_p), step=48))
        ax4.grid(True)
        ax4.set_ylabel(r'累加雨量(mm/39hr)', fontsize=18)
        ax4.legend(loc='lower left', borderaxespad=0, edgecolor='black', facecolor='w', framealpha=1, prop={"family":"MS Gothic"})
        ax4.tick_params(labelsize=18)
        
        fig.text(0.87, 0.905, "【事前放流実施可能性の検討】L2 " + title, ha='right', va='top', color="black", fontsize=25)
        fig.text(0.87, 0.48, 'サーチャージ水位', ha='right', va='top', color="red", fontsize=18)
        fig.text(0.87, 0.459, 'ただし書き操作開始水位', ha='right', va='top', color="orange", fontsize=18)
        fig.text(0.17, 0.937, '———', ha='right', va='top', color="green", fontsize=20)
        fig.text(0.17, 0.917, '－－－', ha='right', va='top', color="green", fontsize=20)
        fig.text(0.17, 0.897, '－・－', ha='right', va='top', color="green", fontsize=20)
        fig.text(0.24, 0.937, '基準降雨量超過', ha='right', va='top', color="black", fontsize=20)
        fig.text(0.24, 0.917, '予備放流開始', ha='right', va='top', color="black", fontsize=20)
        fig.text(0.24, 0.897, '事前放流開始', ha='right', va='top', color="black", fontsize=20)
        plt.rc('legend',fontsize=20)
        plt.tight_layout()
        plt.savefig(f'./tmp.png', dpi=300, bbox_inches="tight", pad_inches=0.05)
        plt.close()


if __name__ == '__main__':
    Hini = 72.0
    extra = 24
    files = glob.glob("./amagase_dam/Qin_L2/*.csv")
    out_list = []
    dam = dam()
    for file in files:
        df_r = read_ame("./amagase_dam/rain_L2/"+os.path.splitext(os.path.basename(file))[0]+".AME", extra)
        shutil.copy(file, "./amagase_dam/Qin.csv")
        dam._reset(Hini,extra)
        dam._loop(df_r, os.path.splitext(os.path.basename(file))[0])
        dam._graph(os.path.splitext(os.path.basename(file))[0])
        shutil.copy("./tmp.png", "./Output_L2/"+os.path.splitext(os.path.basename(file))[0]+".png")
        shutil.copy("./check.csv", "./Output_L2/"+os.path.splitext(os.path.basename(file))[0]+".csv")
