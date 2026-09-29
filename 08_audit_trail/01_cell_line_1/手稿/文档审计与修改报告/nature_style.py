# -*- coding: utf-8 -*-
"""
nature_style.py — 出图统一模块 v2（温和版）
================================================
用法：把本文件与出图脚本放在同一目录，在脚本的 matplotlib 导入之后加一行：

    from nature_style import apply_nature_style
    apply_nature_style()

v2 原则（根据反馈修正）：
  · 不改变你的图布局、不缩尺寸、不删标题——只统一字体与细节样式；
  · 每次 savefig 同时输出两种文件：
      ① .svg（矢量，无限分辨率，供最终出版/排版）
      ② .png（600 dpi，供 Word 插入与屏幕查看——"清晰度不够"就用这个）
  · 初审阶段大图完全合规（Nature 原文："print-publication quality figures are
    large and it is not helpful to upload them at the submission stage"）；
    到最终出版阶段如需 90/180 mm 版心，再调 apply_final_size() 一次即可。
"""
import re
import os
import matplotlib
import matplotlib.pyplot as plt

MM = 1.0 / 25.4
SINGLE_MM, DOUBLE_MM, DEPTH_MM = 90, 180, 170

NATURE_RCPARAMS = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11,                 # 屏幕可读基准；原有大图布局不变
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.titleweight': 'normal',
    'axes.linewidth': 0.8,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'xtick.major.size': 4.0,
    'ytick.major.size': 4.0,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'legend.fontsize': 10,
    'legend.frameon': False,
    'svg.fonttype': 'none',          # SVG 中文字保持可编辑文本
}


def apply_nature_style():
    """应用统一字体与刻度样式（不动布局、不动尺寸）。"""
    matplotlib.rcParams.update(NATURE_RCPARAMS)


def apply_final_size(width_mm=DOUBLE_MM):
    """【最终出版阶段才用】等比收缩到 Nature 版心。初审不要调用。"""
    W, H = width_mm * MM, DEPTH_MM * MM
    for num in plt.get_fignums():
        fig = plt.figure(num)
        w, h = fig.get_size_inches()
        s = min(1.0, W / w, H / h)
        if s < 1.0:
            fig.set_size_inches(w * s, h * s)


# ------------------------------------------------------------------
# savefig 补丁：一次调用，SVG + 600dpi PNG 双输出
# ------------------------------------------------------------------
_orig_savefig = plt.Figure.savefig


def _savefig_dual(self, fname, **kwargs):
    base = re.sub(r'\.(png|jpe?g|tiff?|pdf|svg)$', '', str(fname), flags=re.I)
    kwargs.pop('dpi', None)
    # ① SVG（矢量、文字可编辑）
    _orig_savefig(self, base + '.svg', format='svg', **kwargs)
    # ② PNG（600 dpi，屏幕/Word 用）
    _orig_savefig(self, base + '.png', dpi=600, **kwargs)


plt.Figure.savefig = _savefig_dual
