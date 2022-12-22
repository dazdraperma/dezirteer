#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import time

import gui_support
from sys import platform
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from _version import __version__, __release_year__, __release_month__, __release_date__
import datetime
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as toolbar
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
#from matplotlib.pyplot import errorbar
from matplotlib.collections import LineCollection

from matplotlib.patches import Rectangle
from tkinter import filedialog
from math import sqrt, tan, atan, degrees, cos, sin, pi, floor
from scipy import stats
from numpy import log10, log
from gui_support import set_all_ui_elements

try:
    from Tkinter import *
    from Tkinter import messagebox
except ImportError:
    from tkinter import *
    from tkinter import messagebox
try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True
from math_module import *
from import_export import *
if platform != "darwin":
    import winsound

def truncate(f, n):
    return floor(f * 10 ** n) / 10 ** n

def save_last_dir():
    global g_directory
    dirname=os.getenv('APPDATA')+"\dezirteer"
    completeName = os.path.join(dirname, "lastdir.dzrst")
    f = open(completeName, "w+")
    if isinstance(g_directory, str) and (g_directory != ""):
        f.write(g_directory)   #(g_directory)
    f.close()

def get_last_dir():
    dirname=os.getenv('APPDATA')+"\dezirteer"
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    completeName = os.path.join(dirname, "lastdir.dzrst")
    mode = 'r' if os.path.exists(completeName) else 'a'
    f = open(completeName, mode)
    if mode == "r":
        lines=f.readlines()
        if lines:
            dirpath=lines[0]
        else:
            dirpath = dirname
    else:
        dirpath=dirname
    return dirpath
    f.close()

#calculates KS p- and d-values from current and previous grainsets. If previous grainset doesn't exist,
#sets 0 values to both variables
def set_pval_dval():
    global g_prev_cum, g_prev_prob, g_grainset, g_pval_dval, g_ckde, g_cpdp
    global g_kde, g_pdp
    if g_prev_cum == []:
        pval = 0
        dval = 0

    else:
        if g_graph_settings.pdp_kde_hist == 0:
            curr_cum = g_ckde

        elif g_graph_settings.pdp_kde_hist == 1:
            curr_cum = g_cpdp

        else:
            curr_cum = []
        dval = d_value(curr_cum, g_prev_cum)
        pval = p_value(dval, g_number_of_good_grains[0], g_prev_n[0])

    if g_prev_prob == []:
        like = 0
        sim = 0
    else:
        if g_graph_settings.pdp_kde_hist == 0:
            curr_prob = g_kde
        elif g_graph_settings.pdp_kde_hist == 1:
            curr_prob = g_pdp
        else:
            curr_prob = []
        like = likeness(curr_prob, g_prev_prob)
        sim = similarity(curr_prob, g_prev_prob)
    g_pval_dval = [pval, dval, like, sim]


def peaks():
    global g_kde, g_pdp
    if g_graph_settings.pdp_kde_hist == 0 and len(g_kde) > 0:
        return g_kde[1]

    elif g_graph_settings.pdp_kde_hist == 1 and len(g_pdp) > 0:
        return g_pdp[1]

    else:
        return []


def show_calc_frame(container):
        global g_pval_dval
        frContainer = Frame(container)
        frContainer.configure(relief=GROOVE)
        frContainer.configure(borderwidth="2")
        frContainer.configure(relief=GROOVE)
        frContainer.configure(background="#d9d9d9")
        frContainer.configure(highlightbackground="#d9d9d9")
        frContainer.configure(highlightcolor="black")
        frContainer.pack(fill=BOTH, expand=1)
        container.resizable(False, False)
        elements = ["number of good grains", "weighted average age", "±1σ", "95% conf.", "MSWD", "max age", "min age"]
        list_of_labels = []
        counter = 0
        for n in elements:
            list_of_labels.append(Label(frContainer))
            list_of_labels.append(Label(frContainer))
            list_of_labels[counter*2].grid(row=counter, column=0, pady=5, padx=5, sticky='e')
            list_of_labels[counter*2].configure(text=n)
            list_of_labels[counter*2 + 1].grid(row=counter, column=1, pady=5, padx=5, sticky='w')
            list_of_labels[counter*2 + 1].configure(text=round(g_number_of_good_grains[counter], 2))
            counter += 1

        for x in range(0, 10):
            list_of_labels.append(Label(frContainer))
        list_of_labels[counter * 2 + 0].grid(row=counter, column=0, pady=5, padx=5, sticky='e')
        list_of_labels[counter * 2 + 0].configure(text="peaks: weight")
        list_of_labels[counter * 2 + 1].grid(row=counter, column=1, pady=5, padx=5, sticky='w')
        list_of_labels[counter * 2 + 1].configure(text=calc_peaks_weight(peaks(), g_grainset))

        list_of_labels[counter * 2 + 2].grid(row=counter + 1, column=0, pady=5, padx=5, sticky='e')
        list_of_labels[counter * 2 + 2].configure(text="KS p-val")
        list_of_labels[counter * 2 + 3].grid(row=counter + 1, column=1, pady=5, padx=5, sticky='w')
        list_of_labels[counter * 2 + 3].configure(text=round(g_pval_dval[0], 2))

        list_of_labels[counter * 2 + 4].grid(row=counter + 2, column=0, pady=5, padx=5, sticky='e')
        list_of_labels[counter * 2 + 4].configure(text="KS d-val")
        list_of_labels[counter * 2 + 5].grid(row=counter + 2, column=1, pady=5, padx=5, sticky='w')
        list_of_labels[counter * 2 + 5].configure(text=round(g_pval_dval[1], 2))

        list_of_labels[counter * 2 + 6].grid(row=counter + 3, column=0, pady=5, padx=5, sticky='e')
        list_of_labels[counter * 2 + 6].configure(text="Likeness")
        list_of_labels[counter * 2 + 7].grid(row=counter + 3, column=1, pady=5, padx=5, sticky='w')
        list_of_labels[counter * 2 + 7].configure(text=round(g_pval_dval[2], 2))

        list_of_labels[counter * 2 + 8].grid(row=counter + 4, column=0, pady=5, padx=5, sticky='e')
        list_of_labels[counter * 2 + 8].configure(text="Similarity")
        list_of_labels[counter * 2 + 9].grid(row=counter + 4, column=1, pady=5, padx=5, sticky='w')
        list_of_labels[counter * 2 + 9].configure(text=round(g_pval_dval[3], 2))


class OperationWindow(Frame):
    def __init__(self, master):
        Frame.__init__(self, master)
        self.master = master
        global g_filters, g_grainset, g_list_col_names, g_ckde, g_cpdp, g_kde, g_pdp

        _bgcolor = '#d9d9d9'
        _fgcolor = '#000000'
        _compcolor = '#d9d9d9'
        _ana1color = '#d9d9d9'
        _ana2color = '#d9d9d9'
        font9 = "-family {Segoe UI} -size 8 -weight bold -slant roman" \
                " -underline 0 -overstrike 0"

        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.', background=_bgcolor)
        self.style.configure('.', foreground=_fgcolor)
        self.style.configure('.', font="TkDefaultFont")
        self.style.map('.', background=[('selected', _compcolor), ('active', _ana2color)])

        master.columnconfigure(1, weight=1)
        master.rowconfigure(0, weight=1)
        master.rowconfigure(2, weight=1)

        # _____________________frGraph___________________________________________________________________________________
        self.frGraph = Frame(master)
        self.frGraph.configure(relief=GROOVE)
        self.frGraph.configure(borderwidth="2")
        self.frGraph.configure(relief=GROOVE)
        self.frGraph.configure(background="#d9d9d9")
        self.frGraph.configure(highlightbackground="#d9d9d9")
        self.frGraph.configure(highlightcolor="black")
        self.frGraph.grid(row=0, rowspan=2, columnspan=3, sticky='nswe')
        self.frGraph.columnconfigure(0, weight=1)
        self.frGraph.columnconfigure(1, weight=1)
        self.frGraph.columnconfigure(2, weight=1)
        self.frGraph.rowconfigure(0, weight=1)

        #______________frCon
        self.frConc = Frame(self.frGraph)
        self.frConc.grid(column=0, row=0, sticky='nswe')
        self.frConc.configure(relief=GROOVE)
        self.frConc.configure(borderwidth="2")
        self.frConc.configure(relief=GROOVE)
        self.frConc.configure(background="#d9d9d9")
        self.frConc.configure(highlightbackground="#d9d9d9")
        self.frConc.configure(highlightcolor="black")

        self.fig = Figure(figsize=(4, 2.15), frameon=False)
        # check if MacOs, then plt.figure(), otherwise won't plot
        if platform == "darwin":
            self.fig = plt.figure()
        self.ax_conc = self.fig.add_subplot(111)
        self.ax_conc.axes
        if g_graph_settings.conc_type == 2:
            self.ax_conc.axes.set_yscale("log")

        self.ax_conc.set_xlabel('207Pb/235U')
        self.ax_conc.set_ylabel('206Pb/238U')
        self.ax_conc.set_title('Concordia')
        self.ax_conc.axes.format_coord = lambda x, y: ""
        try:
            self.ax_conc.plot(list(range(0, EarthAge)), graph_to_draw)
        except UnboundLocalError:
            pass

        self.canvas_conc = FigureCanvasTkAgg(self.fig, self.frConc)
        self.canvas_conc.draw()
        self.canvas_conc.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.frConcToolbar = Frame(self.frGraph)
        self.frConcToolbar.grid(column=0, row=1, sticky='ew')
        self.frConcToolbar.configure(relief=GROOVE)
        self.frConcToolbar.configure(borderwidth="2")
        self.frConcToolbar.configure(relief=GROOVE)
        self.frConcToolbar.configure(background="#d9d9d9")
        self.frConcToolbar.configure(highlightbackground="#d9d9d9")
        self.frConcToolbar.configure(highlightcolor="black")
        self.frConcToolbar.configure(width=100)

        # ______________frProb
        self.frProb = Frame(self.frGraph)
        self.frProb.grid(row=0, column=1, sticky='nswe')
        self.frProb.configure(relief=GROOVE)
        self.frProb.configure(borderwidth="2")
        self.frProb.configure(relief=GROOVE)
        self.frProb.configure(background="#d9d9d9")
        self.frProb.configure(highlightbackground="#d9d9d9")
        self.frProb.configure(highlightcolor="black")

        try:
            if (g_graph_settings.pdp_kde_hist == 0) and (g_kde != []):
                graph_to_draw = g_kde

            elif (g_graph_settings.pdp_kde_hist == 1) and (g_pdp != []):
                graph_to_draw = g_pdp[0]
            else:
                pass
        except NameError:
            pass

        self.fig = Figure(figsize=(4, 2.15), frameon=False)
        # check if MacOs, then plt.figure(), otherwise won't plot
        if platform == "darwin":
            self.fig = plt.figure()
        self.ax_prob = self.fig.add_subplot(111)
        self.ax_prob.axes
        self.ax_prob.axes.format_coord = lambda x, y: ""
        self.ax_prob.set_title('KDE/PDP/Histogram')
        self.ax_prob.axes.get_yaxis().set_visible(False)
        try:
            self.ax_prob.plot(list(range(0, EarthAge)), graph_to_draw)
        except UnboundLocalError:
            pass

        self.canvas_prob = FigureCanvasTkAgg(self.fig, self.frProb)
        self.canvas_prob.draw()
        self.canvas_prob.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.frProbToolbar = Frame(self.frGraph)
        self.frProbToolbar.grid(row=1, column=1, sticky='ew')
        self.frProbToolbar.configure(relief=GROOVE)
        self.frProbToolbar.configure(borderwidth="2")
        self.frProbToolbar.configure(relief=GROOVE)
        self.frProbToolbar.configure(background="#d9d9d9")
        self.frProbToolbar.configure(highlightbackground="#d9d9d9")
        self.frProbToolbar.configure(highlightcolor="black")
        self.frProbToolbar.configure(width=100)

        # ______________frCum
        self.frCum = Frame(self.frGraph)
        self.frCum.grid(row=0, column=2, sticky='nswe')
        self.frCum.configure(relief=GROOVE)
        self.frCum.configure(borderwidth="2")
        self.frCum.configure(relief=GROOVE)
        self.frCum.configure(background="#d9d9d9")
        self.frCum.configure(highlightbackground="#d9d9d9")
        self.frCum.configure(highlightcolor="black")
        self.frCum.configure(height=10)

        try:
            if (g_graph_settings.pdp_kde_hist == 0) and (g_ckde != []):
                graph_to_draw = g_ckde

            elif (g_graph_settings.pdp_kde_hist == 1) and (g_cpdp != []):
                graph_to_draw = g_cpdp
            else:
                pass
        except NameError:
            pass

        self.fig = Figure(figsize=(4, 2.15), frameon=False)
        #check if MacOs, then plt.figure(), otherwise won't plot
        if platform == "darwin":
            self.fig = plt.figure()
        self.ax_cum = self.fig.add_subplot(111)
        self.ax_cum.axes.format_coord = lambda x, y: ""
        self.ax_cum.set_title('Cumulative diagrams')
        self.ax_cum.axes.get_yaxis().set_visible(False)

        try:
            self.ax_cum.plot(list(range(0, EarthAge)), graph_to_draw)
        except UnboundLocalError:
            pass

        self.canvas_cum = FigureCanvasTkAgg(self.fig, self.frCum)
        self.canvas_cum.draw()
        self.canvas_cum.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.frCumToolbar = Frame(self.frGraph)
        self.frCumToolbar.grid(row=1, column=2, sticky='ew')
        self.frCumToolbar.configure(relief=GROOVE)
        self.frCumToolbar.configure(borderwidth="2")
        self.frCumToolbar.configure(relief=GROOVE)
        self.frCumToolbar.configure(background="#d9d9d9")
        self.frCumToolbar.configure(highlightbackground="#d9d9d9")
        self.frCumToolbar.configure(highlightcolor="black")
        self.frCumToolbar.configure(width=100)

        #global toolbarConc, toolbarProb, toolbarCum
        toolbarConc = toolbar(self.canvas_conc, self.frConcToolbar)
        toolbarProb = toolbar(self.canvas_prob, self.frProbToolbar)
        toolbarCum = toolbar(self.canvas_cum, self.frCumToolbar)

        # ________frTable_________________________________________________________________________________________________
        self.frTable = Frame(master, height=100)
        self.frTable.configure(relief=GROOVE)
        self.frTable.configure(borderwidth="2")
        self.frTable.configure(relief=GROOVE)
        self.frTable.configure(background="#d9d9d9")
        self.frTable.configure(highlightbackground="#d9d9d9")
        self.frTable.configure(highlightcolor="black")
        self.frTable.grid(row=2, columnspan=3, sticky='nsew')
        self.style.configure('Treeview.Heading', font="TkDefaultFont")

        self.Table = ScrolledTreeView(self.frTable)
        self.Table.place(relx=0.0, rely=0.0, relheight=1.0, relwidth=1.0)
        self.Table['columns'] = g_list_col_names
        self.Table.bind("<Double-1>", self.tableOnDoubleClick)

        # ________frOper_________________________________________________________________________________________________
        self.frOper = Frame(master)
        self.frOper.grid(row=3, columnspan=3, sticky='sew')

        self.frOper.columnconfigure(1, weight=1)

        self.frImport = Frame(self.frOper)
        self.frImport.grid(row=0, column=0, sticky='ns')
        self.frImport.configure(relief=GROOVE)
        self.frImport.configure(borderwidth="2")
        self.frImport.configure(relief=GROOVE)
        self.frImport.configure(background="#d9d9d9")
        self.frImport.configure(highlightbackground="#d9d9d9")
        self.frImport.configure(highlightcolor="black")

        self.lbImport = Label(self.frImport)
        self.lbImport.grid(row=0, columnspan=3, sticky="ew", pady=5)
        self.lbImport.configure(font=font9)
        self.apply_style(self.lbImport)
        self.lbImport.configure(text="Import data")

        self.btnImport = Button(self.frImport, width=14, height=2)
        self.btnImport.grid(row=2, columnspan=3, pady=5)
        self.apply_style(self.btnImport)
        self.btnImport.configure(pady="0")
        self.btnImport.configure(text="Import")
        self.btnImport.configure(command=lambda: self.open_and_load_file())

        self.lbStatus = Label(self.frImport)
        self.lbStatus.grid(row=3, column=0)
        self.apply_style(self.lbStatus)
        self.lbStatus.configure(text='Status:')

        self.lbShowStatus = Label(self.frImport)
        self.lbShowStatus.grid(row=3, column=1, pady=5, padx=5, sticky='w')
        self.apply_style(self.lbShowStatus)
        self.lbShowStatus.configure(relief=SUNKEN)

        self.lbUncType = Label(self.frImport)
        self.lbUncType.grid(row=4, column=0)
        self.apply_style(self.lbUncType)
        self.lbUncType.configure(text='Uncertainty type:')

        self.rbInternal = Radiobutton(self.frImport)
        self.rbInternal.grid(row=4, column=1, sticky='w')
        self.apply_style(self.rbInternal)
        self.rbInternal.configure(font="TkTextFont")
        self.rbInternal.configure(text="Int.")
        self.rbInternal.configure(variable=gui_support.varUncType, value=1)
        self.rbInternal.configure(command=lambda: gui_support.onChange(23, gui_support.varUncType.get(), pars_onChange,
                                                                       self))
        self.rbInternal.select()

        self.rbPropagated = Radiobutton(self.frImport)
        self.rbPropagated.grid(row=5, column=1, sticky='sw')
        self.apply_style(self.rbPropagated)
        self.rbPropagated.configure(text="Prop.")
        self.rbPropagated.configure(variable=gui_support.varUncType, value=2)
        self.rbPropagated.configure(command=lambda: gui_support.onChange(23, gui_support.varUncType.get(),
                                                                         pars_onChange, self))

        '''self.lbSpeedOrPbc1 = Label(self.frImport)
        self.lbSpeedOrPbc1.grid(row=6, column=0)
        self.apply_style(self.lbSpeedOrPbc1)
        self.lbSpeedOrPbc1.configure(text='Calculate Pbc?')

        self.lbSpeedOrPbc2 = Label(self.frImport)
        self.lbSpeedOrPbc2.grid(row=8, columnspan=2, column=0)
        self.apply_style(self.lbSpeedOrPbc2)
        self.lbSpeedOrPbc2.configure(text='Note: Pbc decreases performance')

        self.rbYesSpeedNoPbc = Radiobutton(self.frImport)
        self.rbYesSpeedNoPbc.grid(row=6, column=1, sticky='sw')
        self.apply_style(self.rbYesSpeedNoPbc)
        self.rbYesSpeedNoPbc.configure(text="No Pbc")
        self.rbYesSpeedNoPbc.configure(variable=gui_support.varSpeedOrPbc, value=0)
        self.rbYesSpeedNoPbc.configure(indicatoron=0)
        self.rbYesSpeedNoPbc.configure(command=lambda: gui_support.onChange(31, gui_support.varSpeedOrPbc.get(),
                                                                         pars_onChange, self))

        self.rbNoSpeedYesPbc = Radiobutton(self.frImport)
        self.rbNoSpeedYesPbc.grid(row=7, column=1, sticky='sw')
        self.apply_style(self.rbNoSpeedYesPbc)
        self.rbNoSpeedYesPbc.configure(text="Do Pbc")
        self.rbNoSpeedYesPbc.configure(variable=gui_support.varSpeedOrPbc, value=1)
        self.rbNoSpeedYesPbc.configure(indicatoron=0)
        self.rbNoSpeedYesPbc.configure(command=lambda: gui_support.onChange(31, gui_support.varSpeedOrPbc.get(),
                                                                         pars_onChange, self))'''

        # _______________frSample________________________________________________________________________________________
        self.frSample = Frame(self.frOper)
        self.frSample.grid(row=0, column=1, sticky='nsew')
        self.frSample.configure(relief=GROOVE)
        self.frSample.configure(borderwidth="2")
        self.frSample.configure(relief=GROOVE)
        self.frSample.configure(background="#d9d9d9")
        self.frSample.configure(highlightbackground="#d9d9d9")
        self.frSample.configure(highlightcolor="black")

        self.lbChooseSample = Label(self.frSample)
        self.lbChooseSample.grid(row=0, columnspan=3, sticky="ew", pady=5)
        self.apply_style(self.lbChooseSample)
        self.lbChooseSample.configure(font=font9)
        self.lbChooseSample.configure(text='Choose sample(-s)')

        scrollbar = Scrollbar(self.frSample, orient=VERTICAL)
        self.lboxSamples = Listbox(self.frSample, selectmode='extended', exportselection=0, yscrollcommand=scrollbar.set)
        self.lboxSamples.config(height=10, width=10)
        scrollbar.config(command=self.lboxSamples.yview)
        scrollbar.grid(row=1, column=1, padx=5, sticky="nsw")
        self.lboxSamples.grid(row=1, column=0, sticky="ew", padx=5)
        self.frSample.columnconfigure(0, weight=1)

        # _______________frAgeDisc________________________________________________________________________________________
        self.frAgeDisc = Frame(self.frOper)
        self.frAgeDisc.configure(relief=GROOVE)
        self.frAgeDisc.configure(borderwidth="2")
        self.frAgeDisc.configure(relief=GROOVE)
        self.frAgeDisc.configure(background="#d9d9d9")
        self.frAgeDisc.configure(highlightbackground="#d9d9d9")
        self.frAgeDisc.configure(highlightcolor="black")
        self.frAgeDisc.grid(row=0, column=2, sticky='ns')

        self.lbWhichAge = Label(self.frAgeDisc)
        self.lbWhichAge.grid(row=0, columnspan=3, sticky='ew', pady=5)
        self.apply_style(self.lbWhichAge)
        self.lbWhichAge.configure(font=font9)
        self.lbWhichAge.configure(text='How to calc best age:')
        self.lbWhichAge.configure(state=DISABLED)

        self.entAgeCutoff = Spinbox(self.frAgeDisc, from_=0, to=EarthAge)
        self.entAgeCutoff.grid(row=1, column=1, pady=5, padx=5, sticky='w')
        self.entAgeCutoff.configure(background="white")
        self.entAgeCutoff.configure(disabledforeground="#a3a3a3")
        self.entAgeCutoff.configure(font="TkFixedFont")
        self.entAgeCutoff.configure(foreground="#000000")
        self.entAgeCutoff.configure(insertbackground="black")
        self.entAgeCutoff.configure(textvariable=gui_support.varAgeCutoff)
        self.entAgeCutoff.configure(command=lambda: gui_support.onChange(19, float(self.entAgeCutoff.get()),
                                                                         pars_onChange, self))
        self.entAgeCutoff.bind('<KeyRelease>', (lambda _:gui_support.onChange(19,
                                                                              float(''.join(c for c in self.entAgeCutoff.get() if (c.isdigit() or c =='.'))),
                                                                              pars_onChange, self)))
        self.entAgeCutoff.configure(state=DISABLED)
        self.entAgeCutoff.configure(width=5)

        self.lblAgeMa = Label(self.frAgeDisc)
        self.lblAgeMa.grid(row=1, column=2, sticky='w', pady=5)
        self.apply_style(self.lblAgeMa)
        self.lblAgeMa.configure(text="Ma")

        self.cbWhichAge = ttk.Combobox(self.frAgeDisc)
        self.cbWhichAge.grid(row=1, column=0, sticky='ew')
        self.cbWhichAge.configure(width=15)
        self.cbWhichAge.configure(takefocus="")
        self.cbWhichAge.configure(state=DISABLED)
        self.cbWhichAge.configure(values=('From lesser error', 'Fixed Limit', '207Pb/206Pb', '206Pb/238U', '208Pb/232Th'))
        self.cbWhichAge.current(0)
        self.cbWhichAge.bind('<<ComboboxSelected>>', lambda event: gui_support.onChange(3, self.cbWhichAge.current(),
                                                                                        pars_onChange, self))

        self.lbPbc = Label(self.frAgeDisc)
        self.lbPbc.grid(row=5, sticky='ew', pady=10, columnspan=3)
        self.apply_style(self.lbPbc)
        self.lbPbc.configure(font=font9)
        self.lbPbc.configure(state=DISABLED)
        self.lbPbc.configure(text='Use common Pb corr. ages?')

        self.cbPbc = ttk.Combobox(self.frAgeDisc)
        self.cbPbc.grid(row=6, column=0, sticky='ew')
        self.cbPbc.configure(width=15)
        self.cbPbc.configure(takefocus="")
        self.cbPbc.configure(values=('No', '204Pbc', '207Pbc', '208Pbc', 'Ander.'))
        self.cbPbc.current(0)
        self.cbPbc.bind('<<ComboboxSelected>>', lambda event: gui_support.onChange(4, self.cbPbc.current(),
                                                                                   pars_onChange, self))

        self.entAgeAndersen = Spinbox(self.frAgeDisc, from_=0, to=EarthAge)
        self.entAgeAndersen.grid(row=6, column=1, pady=5, padx=5, sticky='w')
        self.entAgeAndersen.configure(background="white")
        self.entAgeAndersen.configure(disabledforeground="#a3a3a3")
        self.entAgeAndersen.configure(font="TkFixedFont")
        self.entAgeAndersen.configure(foreground="#000000")
        self.entAgeAndersen.configure(insertbackground="black")
        self.entAgeAndersen.configure(textvariable=gui_support.varAgeAndersen)
        self.entAgeAndersen.configure(
            command=lambda: gui_support.onChange(28, float(self.entAgeAndersen.get()), pars_onChange))
        self.entAgeAndersen.bind('<KeyRelease>', (lambda _: gui_support.onChange(28,
                                                                               float(''.join(
                                                                                   c for c in self.entAgeAndersen.get() if
                                                                                   (c.isdigit() or c == '.'))),
                                                                               pars_onChange, self)))
        self.entAgeAndersen.configure(state=DISABLED)
        self.entAgeAndersen.configure(width=5)


        self.lblAgeAndersen = Label(self.frAgeDisc)
        self.lblAgeAndersen.grid(row=6, column=2, sticky='w', pady=5)
        self.apply_style(self.lblAgeAndersen)
        self.lblAgeAndersen.configure(text="And.Ma")


        self.lbCalcDisc = Label(self.frAgeDisc)
        self.lbCalcDisc.grid(row=7, columnspan=2, sticky='ew', pady=5)
        self.apply_style(self.lbCalcDisc)
        self.lbCalcDisc.configure(font=font9)
        self.lbCalcDisc.configure(state=DISABLED)
        self.lbCalcDisc.configure(text='How to calc discordance:')

        self.entDiscAgeFixedLim = Spinbox(self.frAgeDisc, from_=0, to=EarthAge)
        self.entDiscAgeFixedLim.grid(row=8, column=1, pady=5, padx=5, sticky='w')
        self.entDiscAgeFixedLim.configure(background="white")
        self.entDiscAgeFixedLim.configure(disabledforeground="#a3a3a3")
        self.entDiscAgeFixedLim.configure(font="TkFixedFont")
        self.entDiscAgeFixedLim.configure(foreground="#000000")
        self.entDiscAgeFixedLim.configure(insertbackground="black")
        self.entDiscAgeFixedLim.configure(textvariable=gui_support.varDiscCutoff)
        self.entDiscAgeFixedLim.configure(command=lambda: gui_support.onChange(25, float(self.entDiscAgeFixedLim.get()),
                                                                               pars_onChange, self))
        self.entDiscAgeFixedLim.bind('<KeyRelease>', (lambda _: gui_support.onChange(25, float(
            ''.join(c for c in self.entDiscAgeFixedLim.get() if (c.isdigit() or c == '.'))), pars_onChange, self)))
        self.entDiscAgeFixedLim.configure(state=DISABLED)
        self.entDiscAgeFixedLim.configure(width=5)

        self.lblDiscMa = Label(self.frAgeDisc)
        self.lblDiscMa.grid(row=8, column=2, sticky='w', pady=5)
        self.apply_style(self.lblDiscMa)
        self.lblDiscMa.configure(text="Ma")

        self.cbWhichConc = ttk.Combobox(self.frAgeDisc)
        self.cbWhichConc.grid(row=8, column=0, sticky='ew')
        self.cbWhichConc.configure(width=15)
        self.cbWhichConc.configure(takefocus="")
        self.cbWhichConc.configure(state=DISABLED)

        self.cbWhichConc.configure(values=('Fixed limit (Ma):', '207/206-206/238', '207/235-206/238', 'Lesser of 2'))
        self.cbWhichConc.bind('<<ComboboxSelected>>', lambda event: gui_support.onChange(8, self.cbWhichConc.current()+1,
                                                                                        pars_onChange, self))
        self.cbWhichConc.current(3)




        # _______________frFilter_________________________________________________________________________________________
        self.frFilter = Frame(self.frOper)
        self.frFilter.grid(row=0, column=3, sticky='ns')
        self.frFilter.configure(relief=GROOVE)
        self.frFilter.configure(borderwidth="2")
        self.frFilter.configure(relief=GROOVE)
        self.frFilter.configure(background="#d9d9d9")
        self.frFilter.configure(highlightbackground="#d9d9d9")
        self.frFilter.configure(highlightcolor="black")

        self.lbDiscFilt = Label(self.frFilter)
        self.lbDiscFilt.grid(row=0, column=1, columnspan=3, sticky='w', pady=5)
        self.apply_style(self.lbDiscFilt)
        self.lbDiscFilt.configure(font=font9)
        self.lbDiscFilt.configure(text='EITHER: Filter out high discord. values(%)')

        self.rbDiscPerc = Radiobutton(self.frFilter)
        self.rbDiscPerc.grid(row=0, column=0, sticky='e')
        self.apply_style(self.rbDiscPerc)
        self.rbDiscPerc.configure(font="TkTextFont")
        self.rbDiscPerc.configure(variable=gui_support.varDiscPerc, value=0)
        self.rbDiscPerc.configure(command=lambda: gui_support.onChange(29, gui_support.varDiscPerc.get(), pars_onChange,
                                                                       self))
        self.rbDiscPerc.select()

        self.lbDiscIntersect = Label(self.frFilter)
        self.lbDiscIntersect.grid(row=2, column=1, columnspan=4, sticky='w', pady=5)
        self.apply_style(self.lbDiscIntersect)
        self.lbDiscIntersect.configure(font=font9)
        self.lbDiscIntersect.configure(text='OR: Filter out if miss concordia within uncertainty: ')

        self.cbDiscIntersect = ttk.Combobox(self.frFilter)
        self.cbDiscIntersect.grid(row=3, column=1, sticky='w')
        self.cbDiscIntersect.configure(width=3)
        self.cbDiscIntersect.configure(takefocus="")
        self.cbDiscIntersect.configure(state=DISABLED)
        self.cbDiscIntersect.configure(values=('1σ', '2σ', '3σ', '4σ', '5σ', '6σ', '7σ', '8σ', '9σ', '10σ'))
        self.cbDiscIntersect.current(0)
        self.cbDiscIntersect.bind('<<ComboboxSelected>>',
                                  lambda event: gui_support.onChange(30, self.cbDiscIntersect.current(), pars_onChange,
                                                                     self))

        self.rbDiscIntersect = Radiobutton(self.frFilter)
        self.rbDiscIntersect.grid(row=2, column=0, sticky='e')
        self.apply_style(self.rbDiscIntersect)
        self.rbDiscIntersect.configure(font="TkTextFont")
        self.rbDiscIntersect.configure(variable=gui_support.varDiscPerc, value=1)
        self.rbDiscIntersect.configure(command=lambda: gui_support.onChange(29, gui_support.varDiscPerc.get(), pars_onChange,
                                                                       self))


        self.lblMinus = Label(self.frFilter)
        self.lblMinus.grid(row=1, column=0, sticky='ew', pady=5)
        self.apply_style(self.lblMinus)
        self.lblMinus.configure(text="(-)")

        self.lblPlus = Label(self.frFilter)
        self.lblPlus.grid(row=1, column=2, sticky='ew', pady=5)
        self.apply_style(self.lblPlus)
        self.lblPlus.configure(text="(+)")

        self.entNegDiscFilt = Spinbox(self.frFilter, from_=0, to=1000)
        self.entNegDiscFilt.grid(row=1, column=1,  sticky='w')
        self.entNegDiscFilt.configure(background="white")
        self.entNegDiscFilt.configure(disabledforeground="#a3a3a3")
        self.entNegDiscFilt.configure(font="TkFixedFont")
        self.entNegDiscFilt.configure(foreground="#000000")
        self.entNegDiscFilt.configure(insertbackground="black")
        self.entNegDiscFilt.configure(textvariable=gui_support.varNegDiscFilter)
        self.entNegDiscFilt.configure(command=lambda: gui_support.onChange(7, float(self.entNegDiscFilt.get()),
                                                                          pars_onChange, self))
        self.entNegDiscFilt.bind('<KeyRelease>', (lambda _: gui_support.onChange(7, float(
            ''.join(c for c in self.entNegDiscFilt.get() if (c.isdigit() or c == '.'))), pars_onChange, self)))
        self.entNegDiscFilt.configure(state=DISABLED)
        self.entNegDiscFilt.configure(width=3)

        self.entPosDiscFilt = Spinbox(self.frFilter, from_=0, to=1000)
        self.entPosDiscFilt.grid(row=1, column=3, sticky='w')
        self.entPosDiscFilt.configure(background="white")
        self.entPosDiscFilt.configure(disabledforeground="#a3a3a3")
        self.entPosDiscFilt.configure(font="TkFixedFont")
        self.entPosDiscFilt.configure(foreground="#000000")
        self.entPosDiscFilt.configure(insertbackground="black")
        self.entPosDiscFilt.configure(textvariable=gui_support.varPosDiscFilter)
        self.entPosDiscFilt.configure(
            command=lambda: gui_support.onChange(6, float(self.entPosDiscFilt.get()), pars_onChange, self))
        self.entPosDiscFilt.bind('<KeyRelease>', (lambda _: gui_support.onChange(6, float(
            ''.join(c for c in self.entPosDiscFilt.get() if (c.isdigit() or c == '.'))), pars_onChange, self)))
        self.entPosDiscFilt.configure(state=DISABLED)
        self.entPosDiscFilt.configure(width=3)

        self.lbFilterByError = Label(self.frFilter)
        self.lbFilterByError.grid(row=4, column=0, columnspan=2, pady=5, sticky='ew')
        self.apply_style(self.lbFilterByError)
        self.lbFilterByError.configure(font=font9)
        self.lbFilterByError.configure(state=DISABLED)
        self.lbFilterByError.configure(text='Filter by error:')



        self.entErrFilter = Entry(self.frFilter)
        self.entErrFilter.grid(row=5, column=2, padx=5, sticky='w')
        self.entErrFilter.configure(background="white")
        self.entErrFilter.configure(disabledforeground="#a3a3a3")
        self.entErrFilter.configure(font="TkFixedFont")
        self.entErrFilter.configure(foreground="#000000")
        self.entErrFilter.configure(insertbackground="black")
        self.entErrFilter.configure(textvariable=gui_support.varErrFilter)
        self.entErrFilter.bind('<KeyRelease>', (lambda _: gui_support.onChange(20, float(
            ''.join(c for c in self.entErrFilter.get() if (c.isdigit() or c == '.'))), pars_onChange, self)))
        self.entErrFilter.configure(width=3)
        self.entErrFilter.configure(state=DISABLED)

        self.lblErrCutoff = Label(self.frFilter)
        self.lblErrCutoff.grid(row=4, column=2, sticky='w')
        self.apply_style(self.lblErrCutoff)
        self.lblErrCutoff.configure(text="%")

        self.cbErrFilter = ttk.Combobox(self.frFilter)
        self.cbErrFilter.grid(row=5, column=0, columnspan=2, sticky='w')
        self.cbErrFilter.configure(width=15)
        self.cbErrFilter.configure(takefocus="")
        self.cbErrFilter.configure(state=DISABLED)
        self.cbErrFilter.configure(values=('Not used', 'Used'))
        self.cbErrFilter.configure(width=10)
        self.cbErrFilter.current(0)
        self.cbErrFilter.bind('<<ComboboxSelected>>',
                              lambda event: gui_support.onChange(5, self.cbErrFilter.current(), pars_onChange, self))

        self.chbInclude207235Err = Checkbutton(self.frFilter)
        self.chbInclude207235Err.grid(row=6, column=0, columnspan=3, sticky='w', pady=5)
        self.apply_style(self.chbInclude207235Err)
        self.chbInclude207235Err.configure(text="include 7/5 error")
        self.chbInclude207235Err.configure(justify=LEFT)
        self.chbInclude207235Err.configure(state=DISABLED)
        self.chbInclude207235Err.configure(variable=gui_support.varInclude207235Err)
        self.chbInclude207235Err.configure(command=lambda: gui_support.onChange(22,
                                                                                gui_support.varInclude207235Err.get(),
                                                                                pars_onChange, self))

        self.lbUConcFilter = Label(self.frFilter)
        self.lbUConcFilter.grid(row=4, column=3, columnspan=1, pady=4, sticky='w')
        self.apply_style(self.lbUConcFilter)
        self.lbUConcFilter.configure(font=font9)
        self.lbUConcFilter.configure(text="Filter by Uconc:")

        self.entUconcCutoff = Entry(self.frFilter)
        self.entUconcCutoff.grid(row=5, column=4, columnspan=2, pady=5, padx=5, sticky='w')
        self.entUconcCutoff.configure(background="white")
        self.entUconcCutoff.configure(disabledforeground="#a3a3a3")
        self.entUconcCutoff.configure(font="TkFixedFont")
        self.entUconcCutoff.configure(foreground="#000000")
        self.entUconcCutoff.configure(insertbackground="black")
        self.entUconcCutoff.configure(textvariable=gui_support.varUConc)
        self.entUconcCutoff.bind('<KeyRelease>', (lambda _: gui_support.onChange(18, float(
            ''.join(c for c in self.entUconcCutoff.get() if (c.isdigit() or c == '.'))),
                                                                               pars_onChange, self)))
        self.entUconcCutoff.configure(width=5)
        self.entUconcCutoff.configure(state=DISABLED)

        self.lblUConcCutoff = Label(self.frFilter)
        self.lblUConcCutoff.grid(row=4, column=4, sticky='w', pady=5)
        self.apply_style(self.lblUConcCutoff)
        self.lblUConcCutoff.configure(text="ppm")

        self.cbUConc = ttk.Combobox(self.frFilter)
        self.cbUConc.grid(row=5, column=3, sticky='w')
        self.cbUConc.configure(width=15)
        self.cbUConc.configure(takefocus="")
        self.cbUConc.configure(state=DISABLED)
        self.cbUConc.configure(values=('Not used', 'Used'))
        self.cbUConc.bind('<<ComboboxSelected>>', lambda event: gui_support.onChange(2, self.cbUConc.current(),
                                                                                        pars_onChange, self))
        self.cbUConc.configure(width=10)
        self.cbUConc.current(0)



        # _______________frGraphSettings_________________________________________________________________________________
        self.frGraphSettings = Frame(self.frOper)
        self.frGraphSettings.configure(relief=GROOVE)
        self.frGraphSettings.configure(borderwidth="2")
        self.frGraphSettings.configure(relief=GROOVE)
        self.frGraphSettings.configure(background="#d9d9d9")
        self.frGraphSettings.configure(highlightbackground="#d9d9d9")
        self.frGraphSettings.configure(highlightcolor="black")
        self.frGraphSettings.grid(row=0, column=4, sticky="ns")

        self.lbConc = Label(self.frGraphSettings)
        self.lbConc.grid(row=0, columnspan=4, pady=5, sticky='ew')
        self.apply_style(self.lbConc)
        self.lbConc.configure(font=font9)
        self.lbConc.configure(text='Concordia settings:')

        self.lbConcType = Label(self.frGraphSettings)
        self.lbConcType.grid(row=1, column=0, pady=5, sticky='w')
        self.apply_style(self.lbConcType)
        self.lbConcType.configure(text='Conc. type:')

        self.cbConcType = ttk.Combobox(self.frGraphSettings)
        self.cbConcType.grid(row=1, column=1, sticky='w')
        self.cbConcType.configure(width=15)
        self.cbConcType.configure(takefocus="")
        self.cbConcType.configure(state=DISABLED)
        self.cbConcType.configure(values=('Standard', 'Tera-Wass.'))
        self.cbConcType.bind('<<ComboboxSelected>>', lambda event: gui_support.onGraphChange(g_graph_settings, 0,
                                                                                             self.cbConcType.current()))
        self.cbConcType.config(width=8)
        self.cbConcType.current(0)

        self.lbEllipsesAt = Label(self.frGraphSettings)
        self.lbEllipsesAt.grid(row=1, column=2, pady=5, padx=5, sticky='w')
        self.apply_style(self.lbEllipsesAt)
        self.lbEllipsesAt.configure(text='Ellipses at:')

        self.cbEllipsesAt = ttk.Combobox(self.frGraphSettings)
        self.cbEllipsesAt.grid(row=1, column=3, padx=5, sticky='w')
        self.cbEllipsesAt.configure(width=5)
        self.cbEllipsesAt.configure(takefocus="")
        self.cbEllipsesAt.configure(state=DISABLED)
        self.cbEllipsesAt.configure(values=('1σ', '2σ'))
        self.cbEllipsesAt.bind('<<ComboboxSelected>>',
                             lambda event: gui_support.onGraphChange(g_graph_settings, 2, self.cbEllipsesAt.current()+1))
        self.cbEllipsesAt.config(width=3)
        self.cbEllipsesAt.current(0)

        self.lbShowUncorrCorrBothEllipses = Label(self.frGraphSettings)
        self.lbShowUncorrCorrBothEllipses.grid(row=2, column=0, columnspan=1,  sticky='w')
        self.apply_style(self.lbShowUncorrCorrBothEllipses)
        self.lbShowUncorrCorrBothEllipses.configure(text='204Pbc ellipses?')

        self.cbShowUncorrCorrBothEllipses = ttk.Combobox(self.frGraphSettings)
        self.cbShowUncorrCorrBothEllipses.grid(row=2, column=1, columnspan=2, sticky='w', pady=5)
        self.cbShowUncorrCorrBothEllipses.configure(width=15)
        self.cbShowUncorrCorrBothEllipses.configure(takefocus="")
        self.cbShowUncorrCorrBothEllipses.configure(state=DISABLED)
        self.cbShowUncorrCorrBothEllipses.configure(values=('Uncorr.', '204Pb-corr.', 'Both'))
        self.cbShowUncorrCorrBothEllipses.config(width=8)
        self.cbShowUncorrCorrBothEllipses.current(0)

        self.chbIncludeBadEllipses = Checkbutton(self.frGraphSettings)
        self.chbIncludeBadEllipses.grid(row=2, column=2, columnspan=2, sticky='w', pady=5)
        self.apply_style(self.chbIncludeBadEllipses)
        self.chbIncludeBadEllipses.configure(text="show bad spots")
        self.chbIncludeBadEllipses.configure(justify=LEFT)
        self.chbIncludeBadEllipses.configure(state=DISABLED)
        self.chbIncludeBadEllipses.configure(variable=gui_support.varIncludeBadEllipses)

        self.chbShowErrors = Checkbutton(self.frGraphSettings)
        self.chbShowErrors.grid(row=3, column=2, columnspan=2, sticky='w', pady=5)
        self.apply_style(self.chbShowErrors)
        self.chbShowErrors.configure(text="show error bars")
        self.chbShowErrors.configure(justify=LEFT)
        self.chbShowErrors.configure(state=DISABLED)
        self.chbShowErrors.configure(variable=gui_support.varShowErrorBars)


        self.lbDensityPlot = Label(self.frGraphSettings)
        self.lbDensityPlot.grid(row=4, columnspan=3, pady=5, sticky='ew')
        self.apply_style(self.lbDensityPlot)
        self.lbDensityPlot.configure(font=font9)
        self.lbDensityPlot.configure(text='Density plot:')

        self.lbDensityPlotType = Label(self.frGraphSettings)
        self.lbDensityPlotType.grid(row=5, column=0, pady=5, sticky='w')
        self.apply_style(self.lbDensityPlotType)
        self.lbDensityPlotType.configure(text='Type:')

        self.entKDEBandwidth = Spinbox(self.frGraphSettings, from_=1, to=3000)
        self.entKDEBandwidth.grid(row=4, column=3, pady=5, padx=5, sticky='w')
        self.entKDEBandwidth.configure(background="white")
        self.entKDEBandwidth.configure(disabledforeground="#a3a3a3")
        self.entKDEBandwidth.configure(font="TkFixedFont")
        self.entKDEBandwidth.configure(foreground="#000000")
        self.entKDEBandwidth.configure(insertbackground="black")
        self.entKDEBandwidth.configure(textvariable=gui_support.varKDEBandwidth)
        self.entKDEBandwidth.configure(width=5)
        self.entKDEBandwidth.configure(state=DISABLED)
        self.entKDEBandwidth.configure(width=5)
        self.entKDEBandwidth.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 11, float(self.entKDEBandwidth.get())))
        self.entKDEBandwidth.bind('<KeyRelease>', (lambda _: gui_support.onGraphChange(g_graph_settings, 11,
                                                                                       float(''.join(c for c in
                                                                                                     self.entKDEBandwidth.get()
                                                                                                     if (c.isdigit() or c == '.'))))))



        self.entHistBinwidth = Spinbox(self.frGraphSettings, from_=1, to=3000)
        self.entHistBinwidth.grid(row=5, column=3, pady=5, padx=5, sticky='w')
        self.entHistBinwidth.configure(background="white")
        self.entHistBinwidth.configure(disabledforeground="#a3a3a3")
        self.entHistBinwidth.configure(font="TkFixedFont")
        self.entHistBinwidth.configure(foreground="#000000")
        self.entHistBinwidth.configure(insertbackground="black")
        self.entHistBinwidth.configure(textvariable=gui_support.varHistBinwidth)
        self.entHistBinwidth.configure(width=5)
        self.entHistBinwidth.configure(state=DISABLED)
        self.entHistBinwidth.configure(width=5)
        self.entHistBinwidth.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 12, float(self.entHistBinwidth.get())))
        self.entHistBinwidth.bind('<KeyRelease>', (lambda _: gui_support.onGraphChange(g_graph_settings, 12,
                                                                                       float(''.join(c for c in
                                                                                                     self.entHistBinwidth.get()
                                                                                                     if (c.isdigit() or c == '.'))))))


        self.cbDensityPlotType = ttk.Combobox(self.frGraphSettings)
        self.cbDensityPlotType.grid(row=5, column=1, sticky='w')
        self.cbDensityPlotType.configure(takefocus="")
        self.cbDensityPlotType.configure(state=DISABLED)
        self.cbDensityPlotType.configure(values=('KDE', 'PDP', 'Histogram'))
        self.cbDensityPlotType.bind('<<ComboboxSelected>>',
                               lambda event: gui_support.onGraphChange(g_graph_settings, 7,
                                                                       self.cbDensityPlotType.current(),
                                                                       self.entKDEBandwidth, self.entHistBinwidth))
        self.cbDensityPlotType.config(width=7)
        self.cbDensityPlotType.current(0)

        self.lbKDEBandwidth = Label(self.frGraphSettings)
        self.lbKDEBandwidth.grid(row=4, column=2, pady=5, sticky='w')
        self.apply_style(self.lbKDEBandwidth)
        self.lbKDEBandwidth.configure(text='Bandwidth')


        self.lbHistBinwidth = Label(self.frGraphSettings)
        self.lbHistBinwidth.grid(row=5, column=2, pady=5, sticky='w')
        self.apply_style(self.lbHistBinwidth)
        self.lbHistBinwidth.configure(text='Bin width:')

        self.lbAgeCrop = Label(self.frGraphSettings)
        self.lbAgeCrop.grid(row=6, columnspan=4, sticky='ew', pady=5)
        self.apply_style(self.lbAgeCrop)
        self.lbAgeCrop.configure(font=font9)
        self.lbAgeCrop.configure(text='Age crop:')

        self.entAgeMinCrop = Spinbox(self.frGraphSettings, from_=1, to=EarthAge)
        self.entAgeMinCrop.grid(row=7, column=1, pady=5, sticky='w')
        self.entAgeMinCrop.configure(background="white")
        self.entAgeMinCrop.configure(disabledforeground="#a3a3a3")
        self.entAgeMinCrop.configure(font="TkFixedFont")
        self.entAgeMinCrop.configure(foreground="#000000")
        self.entAgeMinCrop.configure(insertbackground="black")
        self.entAgeMinCrop.configure(width=5)

        self.chbMinAgeCrop = Checkbutton(self.frGraphSettings)
        self.chbMinAgeCrop.grid(row=7, column=0, sticky='w', pady=5)
        self.apply_style(self.chbMinAgeCrop)
        self.chbMinAgeCrop.configure(text="Min age:")
        self.chbMinAgeCrop.configure(justify=LEFT)
        self.chbMinAgeCrop.configure(state=DISABLED)
        self.chbMinAgeCrop.configure(variable=gui_support.varMinAgeCrop)
        self.chbMinAgeCrop.configure(command=lambda: gui_support.onChange(26, self.entAgeMinCrop.get(), pars_onChange,
                                                                          self))

        self.entAgeMaxCrop = Spinbox(self.frGraphSettings, from_=1, to=EarthAge)
        self.entAgeMaxCrop.grid(row=7, column=3, pady=5, sticky='w')
        self.entAgeMaxCrop.configure(background="white")
        self.entAgeMaxCrop.configure(disabledforeground="#a3a3a3")
        self.entAgeMaxCrop.configure(font="TkFixedFont")
        self.entAgeMaxCrop.configure(foreground="#000000")
        self.entAgeMaxCrop.configure(insertbackground="black")
        self.entAgeMaxCrop.configure(width=5)

        self.chbMaxAgeCrop = Checkbutton(self.frGraphSettings)
        self.chbMaxAgeCrop.grid(row=7, column=2, sticky='w', pady=5)
        self.apply_style(self.chbMaxAgeCrop)
        self.chbMaxAgeCrop.configure(text="Max age:")
        self.chbMaxAgeCrop.configure(justify=LEFT)
        self.chbMaxAgeCrop.configure(state=DISABLED)
        self.chbMaxAgeCrop.configure(variable=gui_support.varMaxAgeCrop)
        self.chbMaxAgeCrop.configure(command=lambda: gui_support.onChange(27, self.entAgeMaxCrop.get(), pars_onChange,
                                                                          self))



        # _________________frStatus_________________________________________________________________________________________
        self.frStatus = Frame(master)
        self.frStatus.configure(relief=GROOVE)
        self.frStatus.configure(borderwidth="2")
        self.frStatus.configure(relief=GROOVE)
        self.frStatus.configure(background="#d9d9d9")
        self.frStatus.grid(row=4, columnspan=3, sticky='ew')

        self.btnCalcWindow = Button(self.frStatus)
        self.btnCalcWindow.grid(column=4, row=0, rowspan=2, sticky='e', padx=5, pady=6)
        self.apply_style(self.btnCalcWindow)
        self.btnCalcWindow.configure(text="Statistics")
        self.btnCalcWindow.configure(height=2)
        self.btnCalcWindow.configure(width=20)
        self.btnCalcWindow.configure(command=lambda: self.show_frame())

        self.btnClear = Button(self.frStatus)
        self.btnClear.grid(column=3, row=0, rowspan=2, sticky='e', padx=5, pady=6)
        self.apply_style(self.btnClear)
        self.btnClear.configure(text='Clear graph')
        self.btnClear.configure(height=2)
        self.btnClear.configure(width=20)
        self.btnClear.configure(command=lambda: self.clear_graph())

        self.btnExport = Button(self.frStatus)
        self.btnExport.grid(column=2, row=0, rowspan=2, sticky='e', padx=5, pady=6)
        self.apply_style(self.btnExport)
        self.btnExport.configure(text='Export table')
        self.btnExport.configure(width=20)
        self.btnExport.configure(height=2)
        self.btnExport.configure(command=lambda: self.export_dialog())

        self.chbShowCalc = Checkbutton(self.frStatus)
        self.chbShowCalc.grid(row=0, column=5, padx=5, sticky='w')
        self.apply_style(self.chbShowCalc)
        self.chbShowCalc.configure(justify=LEFT)
        self.chbShowCalc.configure(text='Show peaks and stat.')
        self.chbShowCalc.configure(variable=gui_support.varShowCalc)
        self.chbShowCalc.configure(command=lambda: self.plot_text(g_pval_dval[0], g_pval_dval[1]))

        self.chbKeepPrev = Checkbutton(self.frStatus)
        self.chbKeepPrev.grid(row=0, column=6, padx=5, sticky='w')
        self.apply_style(self.chbKeepPrev)
        self.chbKeepPrev.configure(justify=LEFT)
        self.chbKeepPrev.configure(text='Keep prev.')
        self.chbKeepPrev.configure(variable=gui_support.varKeepPrev)
        self.chbKeepPrev.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 13,
                                                                            gui_support.varKeepPrev,
                                                                            self.chbLimitAgeSpectrum))

        self.chbLimitAgeSpectrum = Checkbutton(self.frStatus)
        self.chbLimitAgeSpectrum.grid(row=0, column=7, pady=5, columnspan=2, sticky='w')
        self.apply_style(self.chbLimitAgeSpectrum)
        self.chbLimitAgeSpectrum.configure(justify=LEFT)
        self.chbLimitAgeSpectrum.configure(text='Zoom to ages')
        self.chbLimitAgeSpectrum.configure(variable=gui_support.varLimitAgeSpectrum)
        self.chbLimitAgeSpectrum.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 13,
                                                                                    gui_support.varLimitAgeSpectrum,
                                                                                    self.chbKeepPrev))

        self.btnDraw = Button(self.frStatus)
        self.btnDraw.grid(column=9, row=0, rowspan=2, sticky='e', padx=5, pady=6)
        self.apply_style(self.btnDraw)
        self.btnDraw.configure(text="Plot")
        self.btnDraw.configure(bg='#00ff80')
        self.btnDraw.configure(height=2)
        self.btnDraw.configure(width=20)
        self.btnDraw.configure(command=lambda: self.clear_and_plot())


        #________________Menu___________________________________________________________________________________________
        self.menubar = Menu(master, font="TkMenuFont", bg=_bgcolor, fg=_fgcolor)
        master.configure(menu=self.menubar)
        global mFile
        mFile = Menu(self.menubar, tearoff=False)
        mEdit = Menu(self.menubar, tearoff=False)
        mAbout = Menu(self.menubar, tearoff=False)

        self.menubar.add_cascade(label="File", menu=mFile, underline=0)
#        self.menubar.add_cascade(label="Edit", menu=mEdit, underline=0)
        self.menubar.add_cascade(label="About", menu=mAbout, underline=0)

        #mFile.add_command(label="New Session", underline=0, accelerator="Ctrl+N",
        #                  command=lambda: self.reset_controls(False))

        mFile.add_command(label="Load Session", underline=0, accelerator="Ctrl+O",
                          command=lambda: self.load_session())

        mFile.add_command(label="Save Session", underline=0, accelerator="Ctrl+S",
                          command=lambda: self.save_session())

        mFile.entryconfig(1, state=DISABLED)

        mFile.add_separator()

        #mFile.add_command(label="Export Table", underline=7, accelerator="Ctrl+E+T")
        #mFile.add_command(label="Export Graph", underline=7, accelerator="Ctrl+E+G")
        #mFile.add_separator()
        mFile.add_command(label="Exit", underline=1, command=root.quit, accelerator="Ctrl+Q")
#        mEdit.add_command(label="Undo", accelerator="Ctrl+Z")
#        mEdit.add_command(label="Redo", accelerator="Ctrl+Shift+Z")
#        mEdit.add_separator()
#        mEdit.add_command(label="Settings")

        self.reset_controls(False)

    #_____________Class Methods_________________________________________________________________________________________
    def enable_all_ui_elements(self):

        for var_frame in (self.frImport, self.frAgeDisc, self.frFilter, self.frGraphSettings, self.frStatus):
            for child in var_frame.winfo_children():
                child.configure(state=NORMAL)



    def get_ui_values(self):
        global g_directory
        gui_elements = []
        gui_elements.append(self.lbShowStatus.cget("text")) #0
        gui_elements.append(gui_support.varUncType.get())   #1
        gui_elements.append(self.cbWhichAge.get())          #2
        gui_elements.append(self.entAgeCutoff.get())        #3
        gui_elements.append(self.cbPbc.get())               #4
        gui_elements.append(self.cbWhichConc.get())         #5
        gui_elements.append(self.entDiscAgeFixedLim.get())  #6
        gui_elements.append(self.entNegDiscFilt.get())      #7
        gui_elements.append(self.entPosDiscFilt.get())      #8
        gui_elements.append(self.cbErrFilter.get())         #9
        gui_elements.append(self.entErrFilter.get())        #10
        gui_elements.append(gui_support.varInclude207235Err.get())  #11
        gui_elements.append(self.entUconcCutoff.get())      #12
        gui_elements.append(self.cbUConc.get())             #13
        gui_elements.append(self.cbConcType.get())          #14
        gui_elements.append(self.cbEllipsesAt.get())        #15
        gui_elements.append(self.cbDensityPlotType.get())   #16
        gui_elements.append(self.entHistBinwidth.get())     #17
        gui_elements.append(self.entKDEBandwidth.get())     #18
        gui_elements.append(self.entAgeMinCrop.get())       #19
        gui_elements.append(gui_support.varMinAgeCrop.get())  #20
        gui_elements.append(self.entAgeMaxCrop.get())       #21
        gui_elements.append(gui_support.varShowCalc.get())  #22
        gui_elements.append(gui_support.varKeepPrev.get())  #23
        gui_elements.append(gui_support.varMaxAgeCrop.get())    #24
        gui_elements.append(self.lboxSamples.curselection()) #25
        gui_elements.append(gui_support.varLimitAgeSpectrum.get()) #26
        gui_elements.append(self.lbShowStatus.cget("fg"))    #27
        gui_elements.append(gui_support.varDiscPerc.get())    #28
        gui_elements.append(self.cbDiscIntersect.get())    #29
        gui_elements.append(self.cbShowUncorrCorrBothEllipses.get())  # 30
        gui_elements.append(gui_support.varIncludeBadEllipses.get())  # 31

        # for web version
        gui_elements.append([self.lboxSamples.get(idx) for idx in self.lboxSamples.curselection()])  # 32 select names
        gui_elements.append(self.cbWhichAge.current())  # 33 (2)
        gui_elements.append(self.cbPbc.current())  # 34 (4)
        gui_elements.append(self.cbWhichConc.current())  # 35 (5)
        gui_elements.append(self.cbErrFilter.current())  # 36 (9)
        gui_elements.append(self.cbUConc.current())  # 37 (13)
        gui_elements.append(self.cbConcType.current())  # 38 (14)
        gui_elements.append(self.cbEllipsesAt.current())  # 39 (15)
        gui_elements.append(self.cbDensityPlotType.current())  # 40 (16)
        gui_elements.append(self.cbDiscIntersect.current())  # 41 (29)
        gui_elements.append(self.cbShowUncorrCorrBothEllipses.current())  # 42 (30)

        # again for the desktop version
        gui_elements.append(gui_support.varShowErrorBars.get())  # 43


        return gui_elements



    def set_ui_values(self, args):
        global g_directory
        self.enable_all_ui_elements()

        self.lbShowStatus.configure(text=args[0])
        gui_support.varUncType.set(args[1])
        self.cbWhichAge.set(args[2])

        self.entAgeCutoff.delete(0, END)
        self.entAgeCutoff.insert(0, args[3])

        self.cbPbc.set(args[4])
        self.cbWhichConc.set(args[5])

        self.entDiscAgeFixedLim.delete(0, END)
        self.entDiscAgeFixedLim.insert(0, args[6])

        self.entNegDiscFilt.delete(0, END)
        self.entNegDiscFilt.insert(0, args[7])

        self.entPosDiscFilt.delete(0, END)
        self.entPosDiscFilt.insert(0, args[8])

        self.cbErrFilter.set(args[9])

        self.entErrFilter.delete(0, END)
        self.entErrFilter.insert(0, args[10])

        gui_support.varInclude207235Err.set(args[11])

        self.entUconcCutoff.delete(0, END)
        self.entUconcCutoff.insert(0, args[12])

        self.cbUConc.set(args[13])
        self.cbConcType.set(args[14])
        self.cbEllipsesAt.set(args[15])
        self.cbDensityPlotType.set(args[16])

        self.entHistBinwidth.delete(0, END)
        self.entHistBinwidth.insert(0, args[17])

        self.entKDEBandwidth.delete(0, END)
        self.entKDEBandwidth.insert(0, args[18])

        gui_support.varMinAgeCrop.set(args[20])
        gui_support.varMaxAgeCrop.set(args[24])

        self.entAgeMinCrop.delete(0, END)
        self.entAgeMinCrop.insert(0, args[19])

        self.entAgeMaxCrop.delete(0, END)
        self.entAgeMaxCrop.insert(0, args[21])

        gui_support.varShowCalc.set(args[22])
        gui_support.varKeepPrev.set(args[23])
        self.fill_listbox()

        for index in args[25]:
            self.lboxSamples.selection_set(index)
        gui_support.varLimitAgeSpectrum.set(args[26])
        self.lbShowStatus.configure(fg=args[27])
        gui_support.varDiscPerc.set(args[28])
        self.cbDiscIntersect.set(args[29])
        self.cbShowUncorrCorrBothEllipses.set(args[30])
        gui_support.varIncludeBadEllipses.set(args[31])
        gui_support.varShowErrorBars.set(args[43])


    def show_frame(self):
        winCalc = Toplevel()
        winCalc.resizable(height=None, width=None)
        show_calc_frame(winCalc)

    def apply_style(self, obj):
        obj.configure(activebackground="#f9f9f9")
        obj.configure(activeforeground="black")
        obj.configure(background="#d9d9d9")
        obj.configure(disabledforeground="#a3a3a3")
        obj.configure(foreground="#000000")
        obj.configure(highlightbackground="#d9d9d9")
        obj.configure(highlightcolor="black")

    def open_and_load_file(self, *args):
        global g_directory
        try:
            try:
                global g_plot_txt, g_directory, g_file_type, g_filters, g_list_col_names, g_list_of_samples, \
                    g_grainset, g_number_of_good_grains, pars_onChange
                if g_plot_txt != "":
                    g_plot_txt.remove()
                keep_prev = False
                g_filters.sample_name_filter = []
                use_pbc = gui_support.varMinAgeCrop.get() == 1

                #when run as a main app

                # for unit test
                if args:
                    user_file = args[0]
                    keep_prev = args[1]


                # when module run directly, not imported
                else:
                    user_file = filedialog.askopenfilename(initialdir=g_directory, title="Select file", filetypes=(("Text files", "*.txt"),
                                                                    ("Comma separated values files", "*.csv"),
                                                                    ("All files", "*.*")))

                if user_file != '':
                    if g_grainset != [] and not args:
                        keep_prev = messagebox.askyesno("Keep previous data?", "Keep previous data?")
                    g_directory = os.path.split(user_file)[0]
                    root.title(user_file + ' — Dezirteer: ' + g_dezirteer_version)
                    an_set = []
                    file = imported_file(user_file)
                    g_file_type = file[1]
                    for i in range(1, file[2]):
                        an = file_to_analysis(file, i)
                        an_set.append(an)

                    if keep_prev:
                        an_set = an_set + g_grainset.analyses_list

                    g_grainset = AnalysesSet(an_set, 'set#1')
                    g_grainset.good_bad_sets(g_filters)

                    pars_onChange = [g_filters, self.Table, g_grainset, g_list_col_names]

                    sys.stdout.flush()
                    g_number_of_good_grains = gui_support.fill_data_table(self.Table, g_grainset, g_filters,
                                                                          g_list_col_names)



                    g_list_of_samples = same_sample_set(g_grainset)
                    self.reset_controls(True)
                    #self.set_all_ui_elements()
                    self.clear_prev_or_remove_text()
                else:
                    pass

                save_last_dir()

            except ValueError:
                self.reset_controls(False)
                self.lbShowStatus.configure(state=NORMAL,
                                            text=g_file_type + " data problem\nbetween grains #{}\nand #{}".
                                            format(file_to_analysis(file, i - 1), file_to_analysis(file, i + 1)),
                                            fg="red")

        except FileNotFoundError:
            pass

    #clears the graph after user presses the btClear
    def clear_graph(self):
        global g_plot_txt, g_prev_n, g_prev_cum, g_prev_prob
        g_prev_n = 0
        g_prev_cum = []
        g_prev_prob = []
        self.ax_conc.clear()
        self.ax_prob.clear()
        self.ax_cum.clear()
        self.canvas_conc.draw()
        self.canvas_cum.draw()
        self.canvas_prob.draw()
        self.btnClear.configure(state=DISABLED)
        self.btnCalcWindow.configure(state=DISABLED)
        g_plot_txt = ""

    def export_dialog(self):
        global g_kde, g_ckde, g_pdp, g_cpdp, g_graph_settings
        file_main = filedialog.asksaveasfile(mode='w', defaultextension=".csv", initialdir=g_directory,
                                             filetypes=(("Comma separated values files", "*.csv"),
                                                          ("All files", "*.*")))
        file_prob = os.path.dirname(str(file_main.name)) + '/' + \
                    os.path.splitext(os.path.basename(str(file_main.name)))[0]+'_prob_cum' + '.csv'

        if g_graph_settings.pdp_kde_hist == 0:
            prob = g_kde
            cum = g_ckde
            bandwidth = str(g_graph_settings.bandwidth)
            kde_or_pdp = "KDE"
        elif g_graph_settings.pdp_kde_hist == 1:
            prob = g_pdp
            cum = g_cpdp
            bandwidth = "n/a"
            kde_or_pdp = "PDP"
        else:
            prob = []
            cum = []
            bandwidth = "n/a"
            kde_or_pdp = "Hist"
        gui_support.export_table(g_grainset, g_filters, g_list_col_names, g_graph_settings, file_main, file_prob, prob, cum, bandwidth, kde_or_pdp)

    def save_session(self):

        filename = filedialog.asksaveasfile(mode='w', defaultextension=".dzr", initialdir=g_directory,
                                             filetypes=[("Dezirteer session", "*.dzr")])
        l_ui_values = self.get_ui_values()
        save_object([g_grainset, g_graph_settings, g_filters, l_ui_values], filename.name)


    def fill_listbox(self):
        global g_list_of_samples
        self.lboxSamples.delete(0, END)
        for item in g_list_of_samples:
            self.lboxSamples.insert(END, item.name)

    def load_session(self):
        global g_grainset, g_graph_settings, g_filters, g_list_of_samples, g_number_of_good_grains, g_list_col_names, \
            pars_onChange, mFile
        user_file = filedialog.askopenfilename(initialdir=g_directory, title="Select file",
                                               filetypes=[("Dezirteer session", "*.dzr")])
        loaded_object = load_object(user_file)

        g_grainset = loaded_object[0]
        g_graph_settings = loaded_object[1]
        g_filters = loaded_object[2]
        g_list_of_samples = same_sample_set(g_grainset)

        g_number_of_good_grains = gui_support.fill_data_table(self.Table, g_grainset, g_filters,
                                                              g_list_col_names)

        self.set_ui_values(loaded_object[3])

        set_all_ui_elements(self)
        pars_onChange = [g_filters, self.Table, g_grainset, g_list_col_names]
        mFile.entryconfig(1, state=NORMAL)

    def reset_controls(self, is_data_present):
        global mFile

        if is_data_present:
            mFile.entryconfig(1, state=NORMAL)
            set_all_ui_elements(self)
            self.fill_listbox()
            if self.lboxSamples.get(0) == '':
                status_text = ' data, bad divider'
                status_color = 'red'
            else:
                status_text = '  data OK'
                status_color = 'green'
            self.lbShowStatus.configure(text=g_file_type+status_text, fg=status_color)
        else:
            mFile.entryconfig(1, state=DISABLED)
            self.lboxSamples.delete(0, END)
            for var_frame in (self.frImport, self.frAgeDisc, self.frFilter, self.frGraphSettings, self.frStatus):
                for child in var_frame.winfo_children():
                    child.configure(state=DISABLED)
            self.btnImport.configure(state='normal')
            self.btnCalcWindow.configure(state='disabled')
            self.lbImport.configure(state='normal')
            self.lbShowStatus.configure(text="No Data", fg="red")
            for i in self.Table.get_children():
                self.Table.delete(i)
        global g_plot_txt
        g_plot_txt = ""

    def tableOnDoubleClick(self, event):
        item = self.Table.selection()[0]
        item_name = self.Table.item(item, "text")
        self.clear_and_plot(item_name)

    def min_max_ages(self):
        # choosing age interval based on user's input
        if gui_support.varLimitAgeSpectrum.get() == 1:
            min_age = g_grainset.min_age
            max_age = g_grainset.max_age
            '''min_age = g_number_of_good_grains[6]
            max_age = g_number_of_good_grains[5]'''

            if self.cbConcType.current() == 0:
                min_conc_x = g_grainset.min_207_235
                max_conc_x = g_grainset.max_207_235

                min_conc_y = g_grainset.min_206_238
                max_conc_y = g_grainset.max_206_238
            else:
                min_conc_x = g_grainset.min_238_206
                max_conc_x = g_grainset.max_238_206

                min_conc_y = g_grainset.min_207_206
                max_conc_y = g_grainset.max_207_206

        else:
            min_age = 0
            max_age = EarthAge
            min_conc_x = 0
            min_conc_y = 0
            if self.cbConcType.current() == 0:
                max_conc_x = 100
                max_conc_y = 1.1
            else:
                max_conc_x = 60
                max_conc_y = 0.7
        return [min_age, max_age, min_conc_x, max_conc_x, min_conc_y, max_conc_y]

    #adds or removes text to the cum_plot, depending on the checked state of the chbShowCalc
    def plot_text(self, pval, dval):
        global g_plot_txt, g_conc_txt_green, g_conc_txt_blue, g_conc_txt_black, g_pval_dval

        if gui_support.varShowCalc.get() == 1:
                text_to_show = \
                "n="+str(g_number_of_good_grains[0]) +"\n" \
                "Min age="+str(int(g_number_of_good_grains[6]))+"; "\
                "Max age="+str(int(g_number_of_good_grains[5]))+"\n" \
                "WA age="+str(round((g_number_of_good_grains[1]), 1))+\
                "+-"+str(2 * round((g_number_of_good_grains[2]), 1))+"(2σ int.);\n" \
                "    +-"+str(round((g_number_of_good_grains[3]), 1))+"(95%conf)\n" \
                "MSWD="+str(round(g_number_of_good_grains[4], 2))+"\n" \
                "KS p-value="+str(round(pval, 2))+"; " \
                "d-value="+str(round(dval, 2))+"\n" \
                "Likeness="+str(round(g_pval_dval[2], 2))+"\n" \
                "Similarity=" + str(round(g_pval_dval[3], 2)) + "\n"\
                "Maximum depositional age=" +str(round(g_number_of_good_grains[8][1][0], 1))+"±"+str(round(g_number_of_good_grains[8][1][1], 1))+"\n" \
                "peaks at "
                i = 1
                for p in peaks():
                    if len(peaks()) > 10 and i == 10:
                        text_to_show += "\n (for more peaks click STATISTICS)"
                        break
                    if i < len(peaks()):
                        text_to_show += str(p)+", "
                    else:
                        text_to_show += str(p)
                    if i % 5 == 0 and i < len(peaks()):
                        text_to_show += "\n    "
                    i += 1
        else:
            if g_plot_txt != "":
                g_plot_txt.remove()
            text_to_show = ""

        g_plot_txt = self.ax_cum.text(0.05, 0.10, text_to_show, transform=self.ax_cum.transAxes)
        g_conc_txt_green = self.ax_conc.text(0.65, 0.30, "Green:uncorr.", color="green", transform=self.ax_conc.transAxes)
        g_conc_txt_blue = self.ax_conc.text(0.65, 0.20, "Blue:204Pbc", color="blue", transform=self.ax_conc.transAxes)
        g_conc_txt_black = self.ax_conc.text(0.65, 0.10, "Dotted:bad", color="black", transform=self.ax_conc.transAxes)
        if g_graph_settings.pdp_kde_hist != 2: #if not histogram
            min_max_ages = self.min_max_ages()
            self.plot_peaks(min_max_ages[0], min_max_ages[1])
        self.canvas_cum.draw()
        self.canvas_prob.draw()

    def concordia_type(self):
        # choosing concordia type base on user's input
        if g_graph_settings.conc_type == 0:  # if conventional concordia
            conc_graph_x = [i[1] for i in concordia_table]
            conc_graph_y = [i[0] for i in concordia_table]
            conc_title = "Conventional Concordia"
            conc_graph_xtitle = "207/235"
            conc_graph_ytitle = "206/238"
            xconc = 1
            yconc = 0
        else:  # Tera-Wasserburgh
            conc_graph_x = [(1 / i[0]) for i in concordia_table]
            conc_graph_y = [i[2] for i in concordia_table]
            conc_title = "Tera-Wasserburg Concordia"
            conc_graph_xtitle = "238/206"
            conc_graph_ytitle = "207/206"
            xconc = 3
            yconc = 2
        '''else: #ln T-W
            #for i in concordia_table
            conc_graph_x = [(1 / i[0]) for i in concordia_table]
            #conc_graph_x = [1 / log(i[0]) for i in concordia_table[1:]]
            conc_graph_y = [i[2] for i in concordia_table]
            #conc_graph_y = [log(i[2]) for i in concordia_table[1:]]
            conc_title = "ln Tera-Wasserburg Concordia"
            conc_graph_xtitle = "ln(238/206)"
            conc_graph_ytitle = "ln(207/206)"
            xconc = 3
            yconc = 2'''
        return [conc_graph_x, conc_graph_y, conc_title, conc_graph_xtitle, conc_graph_ytitle, xconc, yconc]

    def kde_pdp_hist(self):
        # choosing kde/pdp/hist based on user input
        global g_ckde, g_cpdp, g_kde, g_pdp, g_prob_graph_to_draw, g_cum_graph_to_draw, g_prob_title, g_cum_title, g_prob_yaxis_title, g_cum_yaxis_title
        if g_graph_settings.pdp_kde_hist == 0:
            g_prob_graph_to_draw = g_kde[0]
            g_cum_graph_to_draw = g_kde[2]
            g_prob_title = "Kernel Density Estimates (KDE)"
            g_cum_title = "Cumulative KDE"
            g_prob_yaxis_title = "relative probability"
            g_cum_yaxis_title = "cumulative probability"
        elif g_graph_settings.pdp_kde_hist == 1:
            g_prob_graph_to_draw = g_pdp[0]
            g_cum_graph_to_draw = g_cpdp
            g_prob_title = "Probability Density Plot (PDP)"
            g_cum_title = "Cumulative PDP"
            g_prob_yaxis_title = "relative probability"
            g_cum_yaxis_title = "cumulative probability"
        else:
            tuple_list = sorted(list(g_grainset.good_set.values()), key=lambda x: x[0])
            g_prob_graph_to_draw = [x[0] for x in tuple_list]
            g_cum_graph_to_draw = []
            g_prob_title = "Histogram"
            g_cum_title = "Cumulative Histogram"
            g_prob_yaxis_title = "# of analyses"
            g_cum_yaxis_title = "# of analyses"
        return[g_prob_graph_to_draw, g_cum_graph_to_draw, g_prob_title, g_cum_title, g_prob_yaxis_title, g_cum_yaxis_title]

    def draw_concordia_ticks(self, xconc, yconc, min_age, max_age):
        if max_age-min_age > 1000:
            step = 500
        elif 500 < max_age-min_age < 1000:
            step = 250
        elif 100 < max_age-min_age < 500:
            step = 50
        elif 50 < max_age-min_age < 100:
            step = 25
        else:
            step = 10
        
        if log10(min_age) >= 2:
            x = -2
        else:
            x = -1
        for t in range(int(truncate(min_age, x)), int(max_age)+step, step):
            if t == 0:
                t += 1
            x = calc_ratio(t)[xconc]
            y = calc_ratio(t)[yconc]
            if g_graph_settings.conc_type == 2:
                x = log(x)
                y = log(y)
            self.ax_conc.plot(x, y, 'ks', markersize=3)
            self.ax_conc.text(x, y, str(t), style='italic')


    def plot_conc_ellipses(self, args):
        # plots ellipses on concordia-discordia diagram
        current_set = [g_grainset.good_set, g_grainset.bad_set]
        which_ellipse_to_plot = self.cbShowUncorrCorrBothEllipses.current()
        if which_ellipse_to_plot == 2:
            j = 2
        else:
            j = 1
        #plot_good_ellipses = gui_support.varIncludeUncorrEllipses.get()
        #plot_204_ellipses = gui_support.varInclude204Ellipses.get()
        plot_bad_ellipses = gui_support.varIncludeBadEllipses.get()
        plot_error_bars = gui_support.varShowErrorBars.get()

        for i in (0, 1):
            for k in range(0, j):
                for zir in current_set[i]:
                    sigma_level = g_graph_settings.ellipses_at
                    if which_ellipse_to_plot == 0 or (which_ellipse_to_plot == 2 and k == 0):
                        corr_coef_75_68 = zir.corr_coef_75_68
                        corr_coef_86_76 = zir.corr_coef_86_76
                        pb207_u235 = zir.pb207_u235
                        pb206_u238 = zir.pb206_u238
                        u238_pb206 = zir.u238_pb206(False)
                        pb207_pb206 = zir.pb207_pb206
                        oval_color = "green"

                    elif which_ellipse_to_plot == 1 or (which_ellipse_to_plot == 2 and k == 1):
                        corr_coef_75_68 = zir.corr_coef_75_68_204
                        corr_coef_86_76 = zir.corr_coef_86_76_204
                        pb207_u235 = zir.rat75_204corr
                        pb206_u238 = zir.rat68_204corr
                        u238_pb206 = zir.u238_pb206(True)
                        pb207_pb206 = zir.rat76_204corr
                        oval_color = "blue"

                    # conventional concordia
                    if g_graph_settings.conc_type == 0:
                        corr_coef = corr_coef_75_68
                        x_conc = pb207_u235[0]  # x-center of the oval
                        y_conc = pb206_u238[0]  # y-center of the oval
                        x_err = pb207_u235[gui_support.varUncType.get()]
                        y_err = pb206_u238[gui_support.varUncType.get()]
                    # Tera-Wasserburg concordia
                    else:
                        corr_coef = corr_coef_86_76
                        x_conc = u238_pb206[0]
                        x_err = u238_pb206[gui_support.varUncType.get()]
                        y_conc = pb207_pb206[0]
                        y_err = pb207_pb206[gui_support.varUncType.get()]
                    '''else: #ln T-W concordia
                        if u238_pb206[0]>0 and pb207_pb206[0]>0:
                            corr_coef = corr_coef_86_76
                            x_conc = log(u238_pb206[0])
                            x_err = 0.434*(u238_pb206[gui_support.varUncType.get()])/(u238_pb206[0])
                            y_conc = (pb207_pb206[0])
                            y_err = 0.434*(pb207_pb206[gui_support.varUncType.get()])/(pb207_pb206[0])
                        else:
                            shall_plot = False'''

                    if (x_conc > 0) and (x_err > 0) and (y_conc > 0) and (y_err > 0):
                        a1 = x_err * corr_coef * sqrt(2) * sigma_level
                        a2 = y_err * corr_coef * sqrt(2) * sigma_level
                        ang = atan(tan(2 * (atan(a2 / a1))) * corr_coef) / 2
                        chi_sq_fact = stats.chi2.ppf(conf_lim(sigma_level), 2)
                        c1 = 2 * (1 - corr_coef ** 2) * chi_sq_fact
                        c2 = 1 / cos(2 * ang)
                        vx = x_err ** 2
                        vy = y_err ** 2
                        test_major_axis = c1 / ((1 + c2) / vx + (1 - c2) / vy)
                        a = sqrt(test_major_axis)
                        test_minor_axis = c1 / ((1 - c2) / vx + (1 + c2) / vy)
                        b = sqrt(test_minor_axis)

                        if i == 1: #bad grains
                            if plot_bad_ellipses == 1 and ((parse_sample_analysis(zir.analysis_name)[0] in g_filters.sample_name_filter) or g_filters.sample_name_filter == []):
                                if args != "":
                                        if zir.analysis_name == args[0]:
                                            oval_fill = True
                                        else:
                                            oval_fill = False
                                #oval_color = 'grey'
                                shall_plot = True
                                line_thickness = 1
                                line_style = ':'
                            else:
                                shall_plot = False

                        else: #good grains
                            if args != "":
                                if zir.analysis_name == args[0]:
                                    oval_fill = True
                                else:
                                    oval_fill = False
                            else:
                                oval_fill = False
                            #oval_color = 'green'
                            shall_plot = True
                            line_thickness = 1
                            line_style = '-'

                        if shall_plot:
                            el = Ellipse(xy=(x_conc, y_conc), width=a * 2, height=b * 2, angle=degrees(ang), color=oval_color,
                                     fill=oval_fill, linewidth=line_thickness, linestyle=line_style)

                            #rect_coords = [x_conc-x_err, y_conc-y_err, x_err*2, y_err*2]

                            self.ax_conc.add_patch(el)
                            if plot_error_bars == 1:

                                self.ax_conc.errorbar(x_conc, y_conc, xerr=x_err*sigma_level, yerr=y_err*sigma_level, xlolims=True, xuplims=True, uplims=True, lolims=True, ecolor='black', fmt='-o')

    def plot_hist(self, min_age, max_age):
        global g_prob_graph_to_draw, g_prob_yaxis_title, g_cum_yaxis_title
        self.ax_prob.axes.get_yaxis().set_visible(True)
        self.ax_cum.axes.get_yaxis().set_visible(True)
        bin_sequence = []
        age = min_age
        bin_width = float(self.entHistBinwidth.get())
        while age < max_age:
            bin_sequence.append(age)
            age += bin_width
        self.ax_prob.hist(g_prob_graph_to_draw, bins=bin_sequence, density=False, cumulative=False)
        g_prob_yaxis_title = g_cum_yaxis_title = "number of grains"
        self.ax_cum.hist(g_prob_graph_to_draw, bins=bin_sequence, density=False, cumulative=True)

    def set_axes(self, conc_title, conc_graph_xtitle, conc_graph_ytitle, conc_graph_x, conc_graph_y, min_age, max_age,
                 min_conc_x, max_conc_x, min_conc_y, max_conc_y):
        # set axis of all graphs
        global g_prob_title, g_cum_title, g_prob_yaxis_title, g_cum_yaxis_title
        self.ax_conc.set_title(conc_title)
        self.ax_conc.set_xlabel(conc_graph_xtitle, labelpad=-16, fontsize=8, position=(0.54, 1e6))
        self.ax_conc.set_ylabel(conc_graph_ytitle, labelpad=-38, fontsize=8)

        self.ax_prob.set_title(g_prob_title)
        self.ax_prob.set_xlabel('Age (Ma)', labelpad=-16, fontsize=8, position=(0.54, 1e6))
        self.ax_prob.set_ylabel(g_prob_yaxis_title)


        self.ax_cum.set_title(g_cum_title)
        self.ax_cum.set_xlabel('Age (Ma)', labelpad=-16, fontsize=8, position=(0.54, 1e6))
        self.ax_cum.set_ylabel(g_cum_yaxis_title)
        #self.ax_cum.set_ylabel(g_cum_yaxis_title, labelpad=-38, fontsize=8)
        self.ax_conc.plot(conc_graph_x, conc_graph_y)

        self.ax_conc.set_xlim(min_conc_x, max_conc_x)
        self.ax_conc.set_ylim(min_conc_y, max_conc_y)
        if g_graph_settings.conc_type == 2:
            self.ax_conc.axes.set_yscale("log")

    def plot_peaks(self, min_age, max_age):
        global g_kde, g_pdp, g_prob_graph_to_draw, g_prob_title
        g_prob_graph_to_draw = self.kde_pdp_hist()[0]

        self.ax_prob.clear()
        self.canvas_prob.draw()
        l_min_age = min_age
        l_max_age = max_age

        if l_min_age % 2 != 0:
            l_min_age -=1

        if l_max_age % 2 != 0:
            l_max_age -=1

        #max_age -= 1
        range_of_ages = range(l_min_age, l_max_age, 2)
        self.ax_prob.plot(list(range_of_ages), g_prob_graph_to_draw[l_min_age//2: l_max_age//2])
        if gui_support.varShowCalc.get() == 1:
            i = 0
            self.ax_prob.set_title(g_prob_title)
            if g_graph_settings.pdp_kde_hist == 0:
                list_peaks = g_kde[1]
            elif g_graph_settings.pdp_kde_hist == 1:
                list_peaks = g_pdp[1]
            else:
                list_peaks = []

            while i < len(list_peaks):
                self.ax_prob.axvline(list_peaks[i], ymin=0, ymax=0.03, color='red')
                i += 1
        self.ax_prob.set_xlabel('Age (Ma)', labelpad=-16, fontsize=8, position=(0.54, 1e6))
        #self.ax_prob.set_ylabel('test', labelpad=-16, fontsize=8, position=(0, 0))
        self.ax_prob.set_title(g_prob_title)

    def prob_cum_plot(self, min_age, max_age):
        global g_prob_graph_to_draw, g_cum_graph_to_draw
        #max_age -= 2
        l_min_age = min_age
        l_max_age = max_age

        if l_min_age % 2 != 0:
            l_min_age -= 1

        if l_max_age % 2 != 0:
            l_max_age -= 1
        range_of_ages = range(l_min_age, l_max_age, 2)
        self.ax_cum.plot(list(range_of_ages), g_cum_graph_to_draw[l_min_age//2: l_max_age//2])
        self.plot_peaks(l_min_age, l_max_age) #ax_prob.plot is done here

    def prob_cum_hist_plot(self, do_hist, min_age, max_age):
        if not do_hist:
            self.ax_prob.axes.get_yaxis().set_visible(False)
            self.ax_cum.axes.get_yaxis().set_visible(False)
            self.prob_cum_plot(min_age, max_age)
        else:
            self.plot_hist(min_age, max_age)

    def clear_prev_or_remove_text(self):
        # clears previous graph, if user chooses to in the chbKeepPrev, else just removes text from cum_plot
        global g_plot_txt
        if gui_support.varKeepPrev.get() == 0:
            self.clear_graph()
        else:
            if g_plot_txt != "":
                g_plot_txt.remove()
        g_plot_txt = ""

    def plot_conc_text_peaks(self,min_age, max_age):
        global g_prev_n, g_prev_cum, g_prev_prob, g_pval_dval, g_ckde, g_cpdp, g_pdp, g_kde

        self.plot_text(g_pval_dval[0], g_pval_dval[1])
        self.canvas_conc.draw()
        #self.canvas_prob.draw() and self.canvas_cum.draw are executed in plot_text
        self.btnClear.configure(state=NORMAL)
        self.btnCalcWindow.configure(state=NORMAL)
        g_prev_n = g_number_of_good_grains
        if g_graph_settings.pdp_kde_hist == 0:
            g_prev_cum = g_ckde
            g_prev_prob = g_kde
        else:
            g_prev_cum = g_cpdp
            g_prev_prob = g_pdp

    def set_plot_types_and_titles(self, kde_pdp_hist):
        global g_prob_graph_to_draw, g_cum_graph_to_draw, g_prob_title, g_cum_title
        g_prob_graph_to_draw = kde_pdp_hist[0]
        g_cum_graph_to_draw = kde_pdp_hist[1]
        g_prob_title = kde_pdp_hist[2]
        g_cum_title = kde_pdp_hist[3]

    #draws the graph based on the data and user settings. Clears the previous graph, or draws on top of it,
    #depending on user settings
    def clear_and_plot(self, *args):
        global g_filters, g_grainset, g_number_of_good_grains, g_plot_txt, g_prev_cum, g_prev_prob, g_prev_n, g_pval_dval
        global g_cpdp, g_ckde, g_kde, g_pdp, g_wa_mda
        g_filters.sample_name_filter = []

        if gui_support.varMinAgeCrop.get() == 1:
            is_editbox_float(self.entAgeMinCrop, '_Filters__minAgeCrop', 0)

        if gui_support.varMaxAgeCrop.get() == 1:
            is_editbox_float(self.entAgeMaxCrop, '_Filters__maxAgeCrop', EarthAge)

        #gets the user-selected items from the listbox

        item_indexes = self.lboxSamples.curselection()
        items = [self.lboxSamples.get(item_indexes) for item_indexes in item_indexes]
        g_filters.sample_name_filter = items
        g_number_of_good_grains = gui_support.fill_data_table(self.Table, g_grainset, g_filters, g_list_col_names)


        #checks if histogram is to be drawn

        do_hist = (g_graph_settings.pdp_kde_hist == 2)

        if g_graph_settings.pdp_kde_hist == 0:
            #start_kde = time.time()
            g_kde = g_grainset.kde(g_graph_settings.bandwidth)
            g_ckde = g_kde[2]
            #end_kde = time.time()
            #total_kde = end_kde - start_kde
            #print("kde: " + str(total_kde))

        elif g_graph_settings.pdp_kde_hist == 1:
            #start_pdp = time.time()
            g_pdp = g_grainset.pdp(gui_support.varUncType.get())
            g_cpdp = g_pdp[2]
            #end_pdp = time.time()
            #total_pdp = end_pdp - start_pdp
            #print("pdp: " + str(total_pdp))

        set_pval_dval()

        # cropping age interval: either full, or cropped from min_age to max_age
        age_lim = self.min_max_ages()
        if gui_support.varMinAgeCrop.get() == 1:
            min_age = int(self.entAgeMinCrop.get())

            min_conc_x = calc_ratio(float(self.entAgeMinCrop.get()))[1]
            min_conc_y = calc_ratio(float(self.entAgeMinCrop.get()))[0]
        else:
            min_age = age_lim[0]
            min_conc_x = age_lim[2]
            min_conc_y = age_lim[4]
        if gui_support.varMaxAgeCrop.get() == 1:
            max_age = int(self.entAgeMaxCrop.get())
            max_conc_x = calc_ratio(float(self.entAgeMaxCrop.get()))[1]
            max_conc_y = calc_ratio(float(self.entAgeMaxCrop.get()))[0]
        else:    
            max_age = age_lim[1]
            max_conc_x = age_lim[3]
            max_conc_y = age_lim[5]

        if min_age < 0:
            min_age = 2
        elif min_age == 0:
            min_age += 2
        elif (min_age > 0) and (min_age % 2 != 0):
            min_age += 1


        self.clear_prev_or_remove_text()

        #FixNeeded. Currently clears the ax_conc every time.
        self.ax_conc.clear()

        #choosing concordia type
        conctype = self.concordia_type()
        conc_graph_x = conctype[0]
        conc_graph_y = conctype[1]
        conc_title = conctype[2]
        conc_graph_xtitle = conctype[3]
        conc_graph_ytitle = conctype[4]
        xconc = conctype[5]
        yconc = conctype[6]

        # choosing kde/pdp/hist

        l_kde_pdp_hist = self.kde_pdp_hist()
        self.set_plot_types_and_titles(l_kde_pdp_hist)


        # set axis of all graphs
        self.set_axes(conc_title, conc_graph_xtitle, conc_graph_ytitle, conc_graph_x, conc_graph_y, min_age, max_age,
                      min_conc_x, max_conc_x, min_conc_y, max_conc_y)

        self.draw_concordia_ticks(xconc, yconc, min_age, max_age)

        if args:
            user_selected_analysis = args
        else:
            user_selected_analysis = ""


        try:
            #plots ellipses on concordia-discordia diagram
            self.plot_conc_ellipses(user_selected_analysis)

            # plotting KDE/CKDE, PDP/CPDP or histogram
            self.prob_cum_hist_plot(do_hist, min_age, max_age)

        #except ValueError:
        #    self.lbShowStatus.configure(text="value error", fg="red")
        #    print ("value error")

        except TypeError:
            self.lbShowStatus.configure(text="type error", fg="red")
            print("type error")

        finally:
            self.plot_conc_text_peaks(min_age, max_age)
            winsound.Beep(2500, 100)




# The following code is added to facilitate the Scrolled widgets
class AutoScroll(object):
    '''Configure the scrollbars for a widget.'''

    def __init__(self, master):
        #  Rozen. Added the try-except clauses so that this class
        #  could be used for scrolled entry widget for which vertical
        #  scrolling is not supported. 5/7/14.
        try:
            vsb = ttk.Scrollbar(master, orient='vertical', command=self.yview)
        except:
            pass
        hsb = ttk.Scrollbar(master, orient='horizontal', command=self.xview)

        #self.configure(yscrollcommand=_autoscroll(vsb),
        #    xscrollcommand=_autoscroll(hsb))
        try:
            self.configure(yscrollcommand=self._autoscroll(vsb))
        except:
            pass
        self.configure(xscrollcommand=self._autoscroll(hsb))

        self.grid(column=0, row=0, sticky='nsew')
        try:
            vsb.grid(column=1, row=0, sticky='ns')
        except:
            pass
        hsb.grid(column=0, row=1, sticky='ew')

        master.grid_columnconfigure(0, weight=1)
        master.grid_rowconfigure(0, weight=1)

        # Copy geometry methods of master  (taken from ScrolledText.py)
        if py3:
            methods = Pack.__dict__.keys() | Grid.__dict__.keys() \
                  | Place.__dict__.keys()
        else:
            methods = Pack.__dict__.keys() + Grid.__dict__.keys() \
                  + Place.__dict__.keys()

        for meth in methods:
            if meth[0] != '_' and meth not in ('config', 'configure'):
                setattr(self, meth, getattr(master, meth))

    @staticmethod
    def _autoscroll(sbar):
        '''Hide and show scrollbar as needed.'''
        def wrapped(first, last):
            first, last = float(first), float(last)
            if first <= 0 and last >= 1:
                sbar.grid_remove()
            else:
                sbar.grid()
            sbar.set(first, last)
        return wrapped

    def __str__(self):
        return str(self.master)


def _create_container(func):
    '''Creates a ttk Frame with a given master, and use this new frame to
    place the scrollbars and the widget.'''
    def wrapped(cls, master, **kw):
        container = ttk.Frame(master)
        return func(cls, container, **kw)
    return wrapped


class ScrolledTreeView(AutoScroll, ttk.Treeview):
    '''A standard ttk Treeview widget with scrollbars that will
    automatically show/hide as needed.'''
    @_create_container
    def __init__(self, master, **kw):
        ttk.Treeview.__init__(self, master, **kw)
        AutoScroll.__init__(self, master)

# checks editboxes for non-numbers; sets g_filters attributes with editbox values
def is_editbox_float(edit_box, to_assign_to, to_replace_with):
    try:
        g_filters.__dict__[to_assign_to] = float(edit_box.get())
    except ValueError:
        g_filters.__dict__[to_assign_to] = to_replace_with
        edit_box.delete(0, END)
        edit_box.insert(0, to_replace_with)

def main():
    global root, g_list_col_names, g_grainset, g_filters, g_graph_settings, prob_fig, prob_subplot
    global g_list_of_samples, g_directory, g_number_of_good_grains, g_prev_cum, g_prev_prob, g_prev_n
    global g_pdp, g_cpdp, g_kde, g_ckde, g_pval_dval, g_dezirteer_version, g_release_date, g_current_date, g_days_since_release
    global g_prob_graph_to_draw, g_cum_graph_to_draw, g_prob_title, g_cum_title, g_prob_yaxis_title, g_cum_yaxis_title
    g_dezirteer_version = __version__
    g_release_date = datetime.date(__release_year__, __release_month__, __release_date__)
    g_current_date = datetime.date.today()
    g_days_since_release = (g_current_date - g_release_date).total_seconds()
    g_days_since_release = int(divmod(g_days_since_release, 86400)[0])
    g_pdp = []
    g_cpdp = []
    g_kde = []
    g_ckde = []
    g_pval_dval = [-1, -1]
    g_prev_cum = []
    g_prev_prob = []
    last_dir = get_last_dir()
    if isinstance(last_dir, str):
        g_directory=last_dir
    else:
        g_directory=os.getenv('APPDATA')+"\dezirteer"
    g_list_col_names = ['232Th/238U', '232/238Err 1s(Int)', '232/238Err 1s(Prop)',
                        '208Pb/232Th', '208/232Err 1s(Int)', '208/232Err 1s(Prop)',
                        '207Pb/206Pb', '207/206Err 1s(Int)', '207/206Err 1s(Prop)',
                        '207Pb/235U', '207/235Err 1s(Int)', '207/235Err 1s(Prop)',
                        '206Pb/238U', '206/238Err 1s(Int)', '206/238Err 1s(Prop)',
                        'corr. coef.75_68', 'corr. coef.86_76',

                        'Uconc (approx. ppm)', 'UconcErr 1s(Int)', 'UconcErr 1s(Prop)',
                        'pbc (approx. ppm)', 'pbcErr 1s(Int)', 'pbcErr 1s(Prop)',
                        '206Pb/204Pb', '206/204Err 1s(Int)', '206/204Err 1s(Prop)',
                        '207Pb/204Pb', '207/204Err 1s(Int)', '207/204Err 1s(Prop)',
                        '208Pb/204Pb', '208/204Err 1s(Int)', '208/204Err 1s(Prop)',
                        '232Th/204Pb', '232/204Err 1s(Int)', '232/204Err 1s(Prop)',
                        '238U/204Pb', '238/204Err 1s(Int)', '238/204Err 1s(Prop)',

                        'Age 208Pb/232Th', 'Age208/232Err 1s(Int)', 'Age208/232Err 1s(Prop)',
                        'Age 207Pb/206Pb', 'Age207/206Err 1s(Int)', 'Age207/206Err 1s(Prop)',
                        'Age 207Pb/235U', 'Age207/235Err 1s(Int)', 'Age207/235Err 1s(Prop)',
                        'Age 206Pb/238U', 'Age206/238Err 1s(Int)', 'Age206/238Err 1s(Prop)',

                        'Corr.type',

                        'Pb204-corr. 68 rat',
                        'Pb204-corr. 68 rat Err 1s(Int)',
                        'Pb204-corr. 68 rat Err 1s(Prop)',

                        'Pb204-corr. 75 rat',
                        'Pb204-corr. 75 rat Err 1s(Int)',
                        'Pb204-corr. 75 rat Err 1s(Prop)',

                        'Pb204-corr. 82 rat',
                        'Pb204-corr. 82 rat Err 1s(Int)',
                        'Pb204-corr. 82 rat Err 1s(Prop)',

                        'Pb204-corr. 76 rat',
                        'Pb204-corr. 76 rat Err 1s(Int)',
                        'Pb204-corr. 76 rat Err 1s(Prop)',

                        'Pb204-corr. 68 age',
                        'Pb204-corr. 68 age Err 1s(Int)',
                        'Pb204-corr. 68 age Err 1s(Prop)',

                        'Pb204-corr. 75 age',
                        'Pb204-corr. 75 age Err 1s(Int)',
                        'Pb204-corr. 75 age Err 1s(Prop)',

                        'Pb204-corr. 82 age',
                        'Pb204-corr. 82 age Err 1s(Int)',
                        'Pb204-corr. 82 age Err 1s(Prop)',

                        'Pb204-corr. 76 age',
                        'Pb204-corr. 76 age Err 1s(Int)',
                        'Pb204-corr. 76 age Err 1s(Prop)',

                        'Pb207-corr. age',
                        'Pb207-corr. age Err 1s(Int)',
                        'Pb207 age corr.Err 1s(Prop)',

                        'Pb208-corr. age',
                        'Pb208-corr. age Err 1s(Int)',
                        'Pb208-corr. age Err 1s(Prop)',

                        'And-corr. age',
                        'And-corr. age Err 1s(Int)',
                        'And-corr. age Err 1s(Prop)',
                        'And. intercept age',


                        'disc. 207/206-206/238', 'disc. 207/235-206/238',
                        'is grain good?', 'best age system',
                        'best age', 'best ageErr 1s',

                        ]
    fill_pbpb_table()
    fill_concordia_table()
    g_filters = Filters()
    g_prev_n = 0
    g_graph_settings = gui_support.GraphSettings()
    root = Tk()
    root.title('Dezirteer: ' + g_dezirteer_version + ', * ' + str(g_days_since_release) + " days ago")
    root.wm_resizable(1, 1)
    gui_support.set_Tk_var()
    master = OperationWindow(root)
    g_grainset = []
    if __name__ == "__main__":
        root.mainloop()


main()