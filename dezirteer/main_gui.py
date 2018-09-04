#addition for checking connection
abc=1

import os
import gui_support
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from tkinter import filedialog
from math import sqrt, tan, atan, degrees, cos, sin, pi
from scipy import stats

try:
    from Tkinter import *
except ImportError:
    from tkinter import *
try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True
from math_module import *

matplotlib.use("TkAgg")


class OperationWindow(Frame):
    def __init__(self, master):
        Frame.__init__(self, master)
        self.master = master
        global g_filters, g_grainset, g_list_col_names

        _bgcolor = '#d9d9d9'
        _fgcolor = '#000000'
        _compcolor = '#d9d9d9'
        _ana1color = '#d9d9d9'
        _ana2color = '#d9d9d9'
        font9 = "-family {Segoe UI} -size 12 -weight bold -slant roman" \
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

        self.frTable = Frame(master, width=955, height=255)
        self.frTable.configure(relief=GROOVE)
        self.frTable.configure(borderwidth="2")
        self.frTable.configure(relief=GROOVE)
        self.frTable.configure(background="#d9d9d9")
        self.frTable.configure(highlightbackground="#d9d9d9")
        self.frTable.configure(highlightcolor="black")
        self.frTable.grid(row=0, column=0, sticky='NWES')
        self.style.configure('Treeview.Heading', font="TkDefaultFont")

        self.Table = ScrolledTreeView(self.frTable)
        self.Table.place(relx=0.0, rely=0.0, relheight=0.96, relwidth=1.0)
        self.Table['columns'] = g_list_col_names
        self.Table.bind("<Double-1>", self.tableOnDoubleClick)


        #pars_onChange = []


        # ________frOper_________________________________________________________________________________________________
        self.frOper = Frame(master)
        self.frOper.grid(row=1, column=0, sticky='sew')

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
        self.lbImport.configure(text="1. Import data")

        self.lbSeparatorType = Label(self.frImport)
        self.lbSeparatorType.grid(row=1, column=0, pady=5, sticky='w')
        self.apply_style(self.lbSeparatorType)
        self.lbSeparatorType.configure(text='Analysis # separator:')

        self.cbSeparatorType = ttk.Combobox(self.frImport)
        self.cbSeparatorType.grid(row=1, column=1, sticky='w')
        self.cbSeparatorType.configure(textvariable=gui_support.varSeparatorType)
        self.cbSeparatorType.configure(width=3)
        self.cbSeparatorType.configure(values=['_', '-', '.', ','], state='readonly')
        self.cbSeparatorType.current(0)

        self.btnImport = Button(self.frImport, width=14, height=2)
        self.btnImport.grid(row=2, columnspan=3, pady=10)
        self.apply_style(self.btnImport)
        self.btnImport.configure(pady="0")
        self.btnImport.configure(text="Import")
        self.btnImport.configure(command=lambda: self.open_and_load_file())

        self.lbStatus = Label(self.frImport)
        self.lbStatus.grid(row=3, column=0)
        self.apply_style(self.lbStatus)
        self.lbStatus.configure(text='''Status:''')

        self.lbShowStatus = Label(self.frImport)
        self.lbShowStatus.grid(row=3, column=1, pady=5, sticky='w')
        self.apply_style(self.lbShowStatus)
        self.lbShowStatus.configure(relief=SUNKEN)

        self.lbUncType = Label(self.frImport)
        self.lbUncType.grid(row=4, column=0)
        self.apply_style(self.lbUncType)
        self.lbUncType.configure(text='''Iolite uncertainty:''')

        self.cbUncType = ttk.Combobox(self.frImport)
        self.cbUncType.grid(row=4, column=1, sticky='w')
        self.cbUncType.configure(textvariable=gui_support.varUncType)
        self.cbUncType.configure(width=10)
        self.cbUncType.configure(values=['Internal', 'Propagated'], state='readonly')
        self.cbUncType.current(0)
        #self.cbUncType.configure(command=lambda: gui_support.onChange(23, gui_support.varUncType.get(), pars_onChange, ))

        self.lbChooseSample = Label(self.frImport)
        self.lbChooseSample.grid(row=5, columnspan=3, sticky="ew", pady=15)
        self.apply_style(self.lbChooseSample)
        self.lbChooseSample.configure(font=font9)
        self.lbChooseSample.configure(text='''2. Choose sample''')

        self.lboxSamples = Listbox(self.frImport, selectmode='extended', exportselection=0)
        self.lboxSamples.grid(row=6, columnspan=3, sticky="ew", padx=5)

        self.lbUConcFilter = Label(self.frImport)
        self.lbUConcFilter.grid(row=7, columnspan=3, pady=4, sticky='ew')
        self.apply_style(self.lbUConcFilter)
        self.lbUConcFilter.configure(font=font9)
        self.lbUConcFilter.configure(text="3. Filter by Uconc?")

        self.rbNoUconc = Radiobutton(self.frImport)
        self.rbNoUconc.grid(row=8, column=0, sticky='w')
        self.apply_style(self.rbNoUconc)
        self.rbNoUconc.configure(font="TkTextFont")
        self.rbNoUconc.configure(text="Don't filter")
        self.rbNoUconc.select()
        self.rbNoUconc.configure(variable=gui_support.varUConc, value=False)
        self.rbNoUconc.configure(command=lambda: gui_support.onChange(2, False, pars_onChange, self.scUconcCutoff))

        self.rbUseUconc = Radiobutton(self.frImport)
        self.rbUseUconc.grid(row=9, column=0, sticky='sw')
        self.apply_style(self.rbUseUconc)
        self.rbUseUconc.configure(text="Cutoff at")
        self.rbUseUconc.configure(variable=gui_support.varUConc, value=True)
        self.rbUseUconc.configure(command=lambda: gui_support.onChange(2, True, pars_onChange, self.scUconcCutoff))

        self.scUconcCutoff = Scale(self.frImport)
        self.scUconcCutoff.grid(row=9, column=1)
        self.scUconcCutoff.configure(activebackground="#d9d9d9")
        self.scUconcCutoff.configure(sliderlength=20)
        self.scUconcCutoff.configure(background="#d9d9d9")
        self.scUconcCutoff.configure(font="TkTextFont")
        self.scUconcCutoff.configure(foreground="#000000")
        self.scUconcCutoff.configure(highlightbackground="#d9d9d9")
        self.scUconcCutoff.configure(highlightcolor="black")
        self.scUconcCutoff.configure(length="100")
        self.scUconcCutoff.configure(orient="horizontal")
        self.scUconcCutoff.configure(tickinterval="0")
        self.scUconcCutoff.configure(to="3000")
        self.scUconcCutoff.configure(bd=2)
        self.scUconcCutoff.set(1000)
        self.scUconcCutoff.configure(state=DISABLED)
        self.scUconcCutoff.configure(troughcolor="#d9d9d9")
        self.scUconcCutoff.configure(command=lambda x: gui_support.onChange(18, self.scUconcCutoff.get(), pars_onChange))

        self.lblPPM = Label(self.frImport)
        self.lblPPM.grid(row=9, column=2, sticky='se', pady=4)
        self.apply_style(self.lblPPM)
        self.lblPPM.configure(text="ppm")

        # _______________frFilter________________________________________________________________________________________
        self.frFilter = Frame(self.frOper)
        self.frFilter.configure(relief=GROOVE)
        self.frFilter.configure(borderwidth="2")
        self.frFilter.configure(relief=GROOVE)
        self.frFilter.configure(background="#d9d9d9")
        self.frFilter.configure(highlightbackground="#d9d9d9")
        self.frFilter.configure(highlightcolor="black")
        self.frFilter.grid(row=0, column=1, sticky='ns')

        self.lbWhichAge = Label(self.frFilter)
        self.lbWhichAge.grid(row=0, columnspan=3, sticky='ew')
        self.apply_style(self.lbWhichAge)
        self.lbWhichAge.configure(font=font9)
        self.lbWhichAge.configure(text='''4. Best age: 7/6 or 6/8?''')
        self.lbWhichAge.configure(state=DISABLED)

        self.rbUBased = Radiobutton(self.frFilter)
        self.rbUBased.configure(variable=gui_support.varAgebased, value=0)
        self.rbUBased.grid(row=1, sticky='w')
        self.apply_style(self.rbUBased)
        self.rbUBased.configure(justify=LEFT)
        self.rbUBased.configure(text='''Based on U-conc.''')
        self.rbUBased.configure(state=DISABLED)
        self.rbUBased.configure(command=lambda: gui_support.onChange(3, 0, pars_onChange, self.rbDiscUbased, self.scAgeCutoff))

        self.rbAgeFixedLim = Radiobutton(self.frFilter)
        self.rbAgeFixedLim.configure(variable=gui_support.varAgebased, value=1)
        self.rbAgeFixedLim.grid(row=2, column=0, sticky='sw', pady=4)
        self.apply_style(self.rbAgeFixedLim)
        self.rbAgeFixedLim.configure(justify=LEFT)
        self.rbAgeFixedLim.configure(text='''Fixed limit (Ma):''')
        self.rbAgeFixedLim.configure(state=DISABLED)
        self.rbAgeFixedLim.configure(command=lambda: gui_support.onChange(3, 1, pars_onChange, self.rbDiscAgeFixedLim,
                                                                          self.scAgeCutoff))
        self.rbAgeFixedLim.select()

        self.rbAge206_207 = Radiobutton(self.frFilter)
        self.rbAge206_207.configure(variable=gui_support.varAgebased, value=2)
        self.rbAge206_207.grid(row=3, sticky='w')
        self.apply_style(self.rbAge206_207)
        self.rbAge206_207.configure(justify=LEFT)
        self.rbAge206_207.configure(text='''206Pb/207Pb''')
        self.rbAge206_207.configure(state=DISABLED)
        self.rbAge206_207.configure(command=lambda: gui_support.onChange(3, 2, pars_onChange, self.rbDisc67_68,
                                                                         self.scAgeCutoff))

        self.rbAge206_238 = Radiobutton(self.frFilter)
        self.rbAge206_238.configure(variable=gui_support.varAgebased, value=3)
        self.rbAge206_238.grid(row=4, sticky='w', pady=10)
        self.apply_style(self.rbAge206_238)
        self.rbAge206_238.configure(justify=LEFT)
        self.rbAge206_238.configure(text='''206Pb/238U''')
        self.rbAge206_238.configure(state=DISABLED)
        self.rbAge206_238.configure(command=lambda: gui_support.onChange(3, 3, pars_onChange, self.rbDisc75_68,
                                                                         self.scAgeCutoff))

        self.scAgeCutoff = Scale(self.frFilter)
        self.scAgeCutoff.grid(row=2, column=1, sticky='ew')
        self.scAgeCutoff.configure(activebackground="#d9d9d9")
        self.scAgeCutoff.configure(sliderlength=20)
        self.scAgeCutoff.configure(background="#d9d9d9")
        self.scAgeCutoff.configure(font="TkTextFont")
        self.scAgeCutoff.configure(foreground="#000000")
        self.scAgeCutoff.configure(highlightbackground="#d9d9d9")
        self.scAgeCutoff.configure(highlightcolor="black")
        self.scAgeCutoff.configure(length="100")
        self.scAgeCutoff.configure(orient="horizontal")
        self.scAgeCutoff.configure(tickinterval="0")
        self.scAgeCutoff.configure(to="3000")
        self.scAgeCutoff.configure(bd=2)
        self.scAgeCutoff.set(1000)
        self.scAgeCutoff.configure(state=DISABLED)
        self.scAgeCutoff.configure(
            command=lambda x: gui_support.onChange(19, self.scAgeCutoff.get(), pars_onChange, self.scDiscAgeFixedLim))
        self.scAgeCutoff.configure(troughcolor="#d9d9d9")

        self.lbPbc = Label(self.frFilter)
        self.lbPbc.grid(row=5, sticky='ew', pady=15, columnspan=3)
        self.apply_style(self.lbPbc)
        self.lbPbc.configure(font=font9)
        self.lbPbc.configure(state=DISABLED)
        self.lbPbc.configure(text='''5. Uncorr. or Pbc?''')

        self.rbUseUncorr = Radiobutton(self.frFilter)
        self.rbUseUncorr.configure(variable=gui_support.varUncorrOrPbc, value=False)
        self.rbUseUncorr.grid(row=6, sticky='w')
        self.apply_style(self.rbUseUncorr)
        self.rbUseUncorr.configure(justify=LEFT)
        self.rbUseUncorr.configure(text='''Uncorr. for Pbc''')
        self.rbUseUncorr.configure(state=DISABLED)
        self.rbUseUncorr.configure(command=lambda: gui_support.onChange(4, False, pars_onChange))
        self.rbUseUncorr.select()

        self.rbUseCorr = Radiobutton(self.frFilter)
        self.rbUseCorr.grid(row=7, sticky='w')
        self.rbUseCorr.configure(variable=gui_support.varUncorrOrPbc, value=True)
        self.apply_style(self.rbUseCorr)
        self.rbUseCorr.configure(justify=LEFT)
        self.rbUseCorr.configure(state=DISABLED)
        self.rbUseCorr.configure(command=lambda: gui_support.onChange(4, True, pars_onChange))
        self.rbUseCorr.configure(text='''Pbc-corr.''')

        self.cbTypePbc = ttk.Combobox(self.frFilter)
        self.cbTypePbc.grid(row=8, sticky='ew')
        self.cbTypePbc.configure(textvariable=gui_support.varTypePbc)
        self.cbTypePbc.configure(width=15)
        self.cbTypePbc.configure(takefocus="")
        self.cbTypePbc.configure(state=DISABLED)

        self.lbFilterByError = Label(self.frFilter)
        self.lbFilterByError.grid(row=9, columnspan=3, pady=15, sticky='ew')
        self.apply_style(self.lbFilterByError)
        self.lbFilterByError.configure(font=font9)
        self.lbFilterByError.configure(state=DISABLED)
        self.lbFilterByError.configure(text='''6. Filter by error?''')

        self.rbNoErrFilter = Radiobutton(self.frFilter)
        self.rbNoErrFilter.configure(variable=gui_support.varErrFilter, value=False)
        self.rbNoErrFilter.grid(row=10, sticky='w')
        self.apply_style(self.rbNoErrFilter)
        self.rbNoErrFilter.configure(justify=LEFT)
        self.rbNoErrFilter.configure(state=DISABLED)
        self.rbNoErrFilter.configure(command=lambda: gui_support.onChange(5, False, pars_onChange,
                                                                          self.chbInclude207235Err, self.scErrFilter))
        self.rbNoErrFilter.configure(text='''Don't filter''')
        self.rbNoErrFilter.select()

        self.rbUseErrFilter = Radiobutton(self.frFilter)
        self.rbUseErrFilter.configure(variable=gui_support.varErrFilter, value=True)
        self.rbUseErrFilter.grid(row=11, column=0, sticky='sw', pady=4)
        self.apply_style(self.rbUseErrFilter)
        self.rbUseErrFilter.configure(justify=LEFT)
        self.rbUseErrFilter.configure(state=DISABLED)
        self.rbUseErrFilter.configure(text='''Cutoff 'best age' at (%): ''')
        self.rbUseErrFilter.configure(command=lambda: gui_support.onChange(5, True, pars_onChange,
                                                                           self.chbInclude207235Err, self.scErrFilter))

        self.chbInclude207235Err = Checkbutton(self.frFilter)
        self.chbInclude207235Err.grid(row=12, column=0, sticky='w', pady=5)
        self.apply_style(self.chbInclude207235Err)
        self.chbInclude207235Err.configure(text="include error in 207/235?")
        self.chbInclude207235Err.configure(justify=LEFT)
        self.chbInclude207235Err.configure(state=DISABLED)
        self.chbInclude207235Err.configure(variable=gui_support.varInclude207235Err)
        self.chbInclude207235Err.configure(command=lambda: gui_support.onChange(22,
                                                                                gui_support.varInclude207235Err.get(),
                                                                                pars_onChange,))

        self.scErrFilter = Scale(self.frFilter)
        self.scErrFilter.grid(row=11, column=1, sticky='ew')
        self.scErrFilter.configure(activebackground="#d9d9d9")
        self.scErrFilter.configure(sliderlength=20)
        self.scErrFilter.configure(background="#d9d9d9")
        self.scErrFilter.configure(font="TkTextFont")
        self.scErrFilter.configure(foreground="#000000")
        self.scErrFilter.configure(highlightbackground="#d9d9d9")
        self.scErrFilter.configure(highlightcolor="black")
        self.scErrFilter.configure(length="100")
        self.scErrFilter.configure(orient="horizontal")
        self.scErrFilter.configure(tickinterval="0")
        self.scErrFilter.configure(from_="1")
        self.scErrFilter.configure(to="50")
        self.scErrFilter.configure(bd=2)
        self.scErrFilter.set(10)
        self.scErrFilter.configure(state=DISABLED)
        self.scErrFilter.configure(troughcolor="#d9d9d9")
        self.scErrFilter.configure(
            command=lambda x: gui_support.onChange(20, self.scErrFilter.get(), pars_onChange, g_list_col_names))

        # _______________frDisc__________________________________________________________________________________________
        self.frDisc = Frame(self.frOper)
        self.frDisc.grid(row=0, column=2, sticky='ns')
        self.frDisc.configure(relief=GROOVE)
        self.frDisc.configure(borderwidth="2")
        self.frDisc.configure(relief=GROOVE)
        self.frDisc.configure(background="#d9d9d9")
        self.frDisc.configure(highlightbackground="#d9d9d9")
        self.frDisc.configure(highlightcolor="black")

        self.lbDiscFilt = Label(self.frDisc)
        self.lbDiscFilt.grid(row=0, columnspan=2, sticky='ew')
        self.apply_style(self.lbDiscFilt)
        self.lbDiscFilt.configure(font=font9)
        self.lbDiscFilt.configure(text='''7. Discord. filters (%)''')

        self.lbPosDiscFilt = Label(self.frDisc)
        self.lbPosDiscFilt.grid(row=1, columnspan=2, sticky='ew')
        self.apply_style(self.lbPosDiscFilt)
        self.lbPosDiscFilt.configure(anchor='w')
        self.lbPosDiscFilt.configure(text='''Positive:''')

        self.scPosDisc = Scale(self.frDisc)
        self.scPosDisc.grid(row=2, columnspan=2, sticky='ew')
        self.scPosDisc.configure(variable=gui_support.varPosDiscFilter)
        self.scPosDisc.configure(activebackground="#d9d9d9")
        self.scPosDisc.configure(background="#d9d9d9")
        self.scPosDisc.configure(font="TkTextFont")
        self.scPosDisc.configure(foreground="#000000")
        self.scPosDisc.configure(highlightbackground="#d9d9d9")
        self.scPosDisc.configure(highlightcolor="black")
        self.scPosDisc.configure(length="50")
        self.scPosDisc.configure(orient="horizontal")
        self.scPosDisc.configure(tickinterval="10.0")
        self.scPosDisc.configure(to="50.0")
        self.scPosDisc.set(20)
        self.scPosDisc.configure(troughcolor="#d9d9d9")
        self.scPosDisc.configure(state=DISABLED)
        self.scPosDisc.configure(command=lambda x: gui_support.onChange(6, self.scPosDisc.get(), pars_onChange))

        self.lbNegDiscFilt = Label(self.frDisc)
        self.lbNegDiscFilt.grid(row=3, columnspan=2, sticky='ew')
        self.apply_style(self.lbNegDiscFilt)
        self.lbNegDiscFilt.configure(anchor='w')
        self.lbNegDiscFilt.configure(state=DISABLED)
        self.lbNegDiscFilt.configure(text='''Negative:''')

        self.scNegDisc = Scale(self.frDisc)
        self.scNegDisc.grid(row=4, columnspan=2, sticky='ew')
        self.scNegDisc.configure(variable=gui_support.varNegDiscFilter)
        self.scNegDisc.configure(activebackground="#d9d9d9")
        self.scNegDisc.configure(background="#d9d9d9")
        self.scNegDisc.configure(font="TkTextFont")
        self.scNegDisc.configure(foreground="#000000")
        self.scNegDisc.configure(from_="-30.0")
        self.scNegDisc.configure(highlightbackground="#d9d9d9")
        self.scNegDisc.configure(highlightcolor="black")
        self.scNegDisc.configure(length="50")
        self.scNegDisc.configure(orient="horizontal")
        self.scNegDisc.configure(tickinterval="10.0")
        self.scNegDisc.configure(to="0.0")
        self.scNegDisc.set(-10)
        self.scNegDisc.configure(state=DISABLED)
        self.scNegDisc.configure(troughcolor="#d9d9d9")
        self.scNegDisc.configure(command=lambda x: gui_support.onChange(7, self.scNegDisc.get(), pars_onChange))

        self.lbCalcDisc = Label(self.frDisc)
        self.lbCalcDisc.grid(row=5, columnspan=2, sticky='ew', pady=5)
        self.apply_style(self.lbCalcDisc)
        self.lbCalcDisc.configure(font=font9)
        self.scNegDisc.configure(state=DISABLED)
        self.lbCalcDisc.configure(text='''8. Discord. type''')

        self.chbDiscLinked2Age = Checkbutton(self.frDisc)
        self.chbDiscLinked2Age.grid(row=6, columnspan=2, sticky='w', pady=5)
        self.apply_style(self.chbDiscLinked2Age)
        self.chbDiscLinked2Age.configure(text="Linked to the choice in #4")
        self.chbDiscLinked2Age.configure(justify=LEFT)
        self.chbDiscLinked2Age.configure(state=DISABLED)
        self.chbDiscLinked2Age.configure(variable=gui_support.varDiscLinked2Age)
        self.chbDiscLinked2Age.configure(command=lambda: gui_support.onChange(21, gui_support.varDiscLinked2Age.get(),
                                                                              pars_onChange,
                                                                              self.rbDiscUbased,
                                                                              self.rbDiscAgeFixedLim,
                                                                              self.rbDisc67_68,
                                                                              self.rbDisc75_68,
                                                                              self.scDiscAgeFixedLim,
                                                                              self.scAgeCutoff))

        self.rbDiscUbased = Radiobutton(self.frDisc)
        self.rbDiscUbased.configure(variable=gui_support.varDiscType, value=0)
        self.rbDiscUbased.grid(row=7, sticky='w', pady=5)
        self.apply_style(self.rbDiscUbased)
        self.rbDiscUbased.configure(justify=LEFT)
        self.rbDiscUbased.configure(text='''Based on U-conc''')
        self.rbDiscUbased.configure(state=DISABLED)
        self.rbDiscUbased.configure(command=lambda: gui_support.onChange(8, 0, pars_onChange))

        self.rbDiscAgeFixedLim = Radiobutton(self.frDisc)
        self.rbDiscAgeFixedLim.configure(variable=gui_support.varDiscType, value=1)
        self.rbDiscAgeFixedLim.grid(row=8, column=0, sticky='ws', pady=5)
        self.apply_style(self.rbDiscAgeFixedLim)
        self.rbDiscAgeFixedLim.configure(justify=LEFT)
        self.rbDiscAgeFixedLim.configure(text='''Fixed limit (Ma):''')
        self.rbDiscAgeFixedLim.configure(state=DISABLED)
        self.rbDiscAgeFixedLim.configure(command=lambda: gui_support.onChange(8, 1, pars_onChange))
        self.rbDiscAgeFixedLim.select()

        self.scDiscAgeFixedLim = Scale(self.frDisc)
        self.scDiscAgeFixedLim.grid(row=8, column=1, sticky='w')
        self.scDiscAgeFixedLim.configure(activebackground="#d9d9d9")
        self.scDiscAgeFixedLim.configure(sliderlength=20)
        self.scDiscAgeFixedLim.configure(background="#d9d9d9")
        self.scDiscAgeFixedLim.configure(font="TkTextFont")
        self.scDiscAgeFixedLim.configure(foreground="#000000")
        self.scDiscAgeFixedLim.configure(highlightbackground="#d9d9d9")
        self.scDiscAgeFixedLim.configure(highlightcolor="black")
        self.scDiscAgeFixedLim.configure(length="100")
        self.scDiscAgeFixedLim.configure(orient="horizontal")
        self.scDiscAgeFixedLim.configure(tickinterval="0")
        self.scDiscAgeFixedLim.configure(to="3000")
        self.scDiscAgeFixedLim.configure(bd=2)
        self.scDiscAgeFixedLim.set(1000)
        self.scDiscAgeFixedLim.configure(state=DISABLED)
        self.scDiscAgeFixedLim.configure(command=lambda x: gui_support.onChange(19, self.scDiscAgeFixedLim.get(),
                                                                                pars_onChange, self.scDiscAgeFixedLim))
        self.scDiscAgeFixedLim.configure(troughcolor="#d9d9d9")

        self.rbDisc67_68 = Radiobutton(self.frDisc)
        self.rbDisc67_68.configure(variable=gui_support.varDiscType, value=2)
        self.rbDisc67_68.grid(row=9, sticky='w', pady=5)
        self.apply_style(self.rbDisc67_68)
        self.rbDisc67_68.configure(justify=LEFT)
        self.rbDisc67_68.configure(text='''206/207-238/206''')
        self.rbDisc67_68.configure(state=DISABLED)
        self.rbDisc67_68.configure(command=lambda: gui_support.onChange(8, 2, pars_onChange))

        self.rbDisc75_68 = Radiobutton(self.frDisc)
        self.rbDisc75_68.configure(variable=gui_support.varDiscType, value=3)
        self.rbDisc75_68.grid(row=10, sticky='sw', pady=5)
        self.apply_style(self.rbDisc75_68)
        self.rbDisc75_68.configure(justify=LEFT)
        self.rbDisc75_68.configure(text='''235/207-238/206''')
        self.rbDisc75_68.configure(state=DISABLED)
        self.rbDisc75_68.configure(command=lambda: gui_support.onChange(8, 3, pars_onChange))

        # _______________frGraphSettings_________________________________________________________________________________
        self.frGraphSettings = Frame(self.frOper)
        self.frGraphSettings.configure(relief=GROOVE)
        self.frGraphSettings.configure(borderwidth="2")
        self.frGraphSettings.configure(relief=GROOVE)
        self.frGraphSettings.configure(background="#d9d9d9")
        self.frGraphSettings.configure(highlightbackground="#d9d9d9")
        self.frGraphSettings.configure(highlightcolor="black")
        self.frGraphSettings.grid(row=0, column=3, sticky="ns")

        self.lbConc = Label(self.frGraphSettings)
        self.lbConc.grid(row=0, columnspan=4, sticky='ew')
        self.apply_style(self.lbConc)
        self.lbConc.configure(font=font9)
        self.lbConc.configure(text='''9. Concordia''')

        self.lbConcType = Label(self.frGraphSettings)
        self.lbConcType.grid(row=1, column=0, pady=5, sticky='w')
        self.apply_style(self.lbConcType)
        self.lbConcType.configure(text='''Conc.type:''')

        self.rbStdConc = Radiobutton(self.frGraphSettings)
        self.rbStdConc.configure(variable=gui_support.varConcType, value=0)
        self.rbStdConc.grid(row=1, column=1, pady=5, sticky='w')
        self.apply_style(self.rbStdConc)
        self.rbStdConc.configure(justify=LEFT)
        self.rbStdConc.configure(text='''Standard''')
        self.rbStdConc.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 0, 0))
        self.rbStdConc.select()

        self.rbTerWassConc = Radiobutton(self.frGraphSettings)
        self.rbTerWassConc.configure(variable=gui_support.varConcType, value=1)
        self.rbTerWassConc.grid(row=1, column=2, pady=5, sticky='w')
        self.apply_style(self.rbTerWassConc)
        self.rbTerWassConc.configure(justify=LEFT)
        self.rbTerWassConc.configure(text='''Ter-Wass.''')
        self.rbTerWassConc.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 0, 1))

        self.lbEclipsesAt = Label(self.frGraphSettings)
        self.lbEclipsesAt.grid(row=3, column=0, pady=5, sticky='w')
        self.apply_style(self.lbEclipsesAt)
        self.lbEclipsesAt.configure(text='''Eclipses at:''')

        self.rbEcl1Sigma = Radiobutton(self.frGraphSettings)
        self.rbEcl1Sigma.configure(variable=gui_support.varEclipseSigma, value=1)
        self.rbEcl1Sigma.grid(row=3, column=1, columnspan=4, pady=5, sticky='w')
        self.apply_style(self.rbEcl1Sigma)
        self.rbEcl1Sigma.configure(justify=LEFT)
        self.rbEcl1Sigma.configure(text='''1σ''')
        self.rbEcl1Sigma.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 2, 1))
        self.rbEcl1Sigma.select()

        self.rbEcl2Sigma = Radiobutton(self.frGraphSettings)
        self.rbEcl2Sigma.configure(variable=gui_support.varEclipseSigma, value=2)
        self.rbEcl2Sigma.grid(row=3, column=2, pady=5, sticky='w')
        self.apply_style(self.rbEcl2Sigma)
        self.rbEcl2Sigma.configure(justify=LEFT)
        self.rbEcl2Sigma.configure(text='''2σ''')
        self.rbEcl2Sigma.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 2, 2))

        self.chbFitDiscordia = Checkbutton(self.frGraphSettings)
        self.chbFitDiscordia.grid(row=5, pady=5, sticky='w')
        self.apply_style(self.chbFitDiscordia)
        self.chbFitDiscordia.configure(justify=LEFT)
        self.chbFitDiscordia.configure(text="Fit discordia?")
        self.chbFitDiscordia.configure(variable=gui_support.varFitDiscordia)
        self.chbFitDiscordia.configure(command=lambda:
        gui_support.onChange(11, gui_support.varFitDiscordia.get(), pars_onChange))

        self.chbAnchored = Checkbutton(self.frGraphSettings)
        self.chbAnchored.grid(row=6, column=0, pady=5, sticky='w')
        self.apply_style(self.chbAnchored)
        self.chbAnchored.configure(justify=LEFT)
        self.chbAnchored.configure(text='Anchored(Ma):')
        self.chbAnchored.configure(variable=gui_support.varAnchored)

        self.entAnchoredAge = Entry(self.frGraphSettings)
        self.entAnchoredAge.grid(row=6, column=1, columnspan=2, pady=5, sticky='w')
        self.entAnchoredAge.configure(background="white")
        self.entAnchoredAge.configure(disabledforeground="#a3a3a3")
        self.entAnchoredAge.configure(font="TkFixedFont")
        self.entAnchoredAge.configure(foreground="#000000")
        self.entAnchoredAge.configure(insertbackground="black")
        self.entAnchoredAge.configure(width=5)

        self.lbKdePdpHist = Label(self.frGraphSettings)
        self.lbKdePdpHist.grid(row=7, columnspan=3, pady=15, sticky='ew')
        self.apply_style(self.lbKdePdpHist)
        self.lbKdePdpHist.configure(font=font9)
        self.lbKdePdpHist.configure(text='''10. KDE/PDP/Hist''')

        self.rbDrawKDE = Radiobutton(self.frGraphSettings)
        self.rbDrawKDE.configure(variable=gui_support.var_pdp_kde_hist, value=0)
        self.rbDrawKDE.grid(row=8, column=0, pady=5, sticky='e')
        self.apply_style(self.rbDrawKDE)
        self.rbDrawKDE.configure(text='''KDE''')
        self.rbDrawKDE.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 7, 0))
        self.rbDrawKDE.select()

        self.rbDrawPdp = Radiobutton(self.frGraphSettings)
        self.rbDrawPdp.configure(variable=gui_support.var_pdp_kde_hist, value=1)
        self.rbDrawPdp.grid(row=8, column=1, pady=5, sticky='ew')
        self.apply_style(self.rbDrawPdp)
        self.rbDrawPdp.configure(text='''PDP''')
        self.rbDrawPdp.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 8, 1))

        self.rbDrawHist = Radiobutton(self.frGraphSettings)
        self.rbDrawHist.configure(variable=gui_support.var_pdp_kde_hist, value=2)
        self.rbDrawHist.grid(row=8, column=2, pady=5, sticky='w')
        self.apply_style(self.rbDrawHist)
        self.rbDrawHist.configure(text='''Hist.''')
        self.rbDrawHist.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 8, 2))

        self.scBandwidth = Scale(self.frGraphSettings)
        self.scBandwidth.grid(row=9, columnspan=4, sticky='ew')
        self.scBandwidth.configure(activebackground="#d9d9d9")
        self.scBandwidth.configure(background="#d9d9d9")
        self.scBandwidth.configure(font="TkTextFont")
        self.scBandwidth.configure(foreground="#000000")
        self.scBandwidth.configure(highlightbackground="#d9d9d9")
        self.scBandwidth.configure(highlightcolor="black")
        self.scBandwidth.configure(length="20")
        self.scBandwidth.configure(label="KDE's bandwidth (Ma):")
        self.scBandwidth.configure(orient="horizontal")
        self.scBandwidth.configure(from_="1")
        self.scBandwidth.configure(tickinterval="499")
        self.scBandwidth.configure(to="500")
        self.scBandwidth.configure(bd=2)
        self.scBandwidth.set(50)
        self.scBandwidth.configure(troughcolor="#d9d9d9")
        self.scBandwidth.configure(command=lambda x: gui_support.onGraphChange(g_graph_settings, 11,
                                                                                       self.scBandwidth.get()))

        self.scBinwidth = Scale(self.frGraphSettings)
        self.scBinwidth.grid(row=10, columnspan=4, sticky='ew')
        self.scBinwidth.configure(activebackground="#d9d9d9")
        self.scBinwidth.configure(background="#d9d9d9")
        self.scBinwidth.configure(font="TkTextFont")
        self.scBinwidth.configure(foreground="#000000")
        self.scBinwidth.configure(highlightbackground="#d9d9d9")
        self.scBinwidth.configure(highlightcolor="black")
        self.scBinwidth.configure(length="20")
        self.scBinwidth.configure(label="Histogram bin width (Ma):")
        self.scBinwidth.configure(orient="horizontal")
        self.scBinwidth.configure(from_="1")
        self.scBinwidth.configure(tickinterval="999")
        self.scBinwidth.configure(to="1000")
        self.scBinwidth.configure(bd=2)
        self.scBinwidth.set(50)
        self.scBinwidth.configure(troughcolor="#d9d9d9")
        self.scBinwidth.configure(command=lambda x: gui_support.onGraphChange(g_graph_settings, 12,
                                                                                      self.scBinwidth.get()))

        self.cbKeepPrev = Checkbutton(self.frGraphSettings)
        self.cbKeepPrev.grid(row=11, column=0, pady=5, sticky='ew')
        self.apply_style(self.cbKeepPrev)
        self.cbKeepPrev.configure(justify=LEFT)
        self.cbKeepPrev.configure(text='''Keep prev.''')
        self.cbKeepPrev.configure(variable=gui_support.varKeepPrev)
        self.cbKeepPrev.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 13,
                                                                              gui_support.varKeepPrev,
                                                                              self.cbLimitAgeSpectrum))

        self.cbLimitAgeSpectrum = Checkbutton(self.frGraphSettings)
        self.cbLimitAgeSpectrum.grid(row=11, column=1, pady=5, sticky='ew')
        self.apply_style(self.cbLimitAgeSpectrum)
        self.cbLimitAgeSpectrum.configure(justify=LEFT)
        self.cbLimitAgeSpectrum.configure(text='''Zoom to ages''')
        self.cbLimitAgeSpectrum.configure(variable=gui_support.varLimitAgeSpectrum)
        self.cbLimitAgeSpectrum.configure(command=lambda: gui_support.onGraphChange(g_graph_settings, 13,
                                                                              gui_support.varLimitAgeSpectrum,
                                                                              self.cbKeepPrev))

        # _________________frStatus_________________________________________________________________________________________
        self.frStatus = Frame(master)
        self.frStatus.configure(relief=GROOVE)
        self.frStatus.configure(borderwidth="2")
        self.frStatus.configure(relief=GROOVE)
        self.frStatus.configure(background="#d9d9d9")
        self.frStatus.grid(row=2, column=0, sticky='ew')

        self.btnDraw = Button(self.frStatus)
        self.btnDraw.grid(column=4, row=0, sticky='e', padx=5, pady=6)
        self.apply_style(self.btnDraw)
        self.btnDraw.configure(text="(Re-)Draw")
        self.btnDraw.configure(height=2)
        self.btnDraw.configure(width=20)
        self.btnDraw.configure(command=lambda: self.clear_and_plot())

        self.btnClear = Button(self.frStatus)
        self.btnClear.grid(column=3, row=0, sticky='e', padx=5, pady=6)
        self.apply_style(self.btnClear)
        self.btnClear.configure(text='''Clear graph''')
        self.btnClear.configure(height=2)
        self.btnClear.configure(width=20)
        self.btnClear.configure(command=lambda: self.clear_graph())

        self.btnExport = Button(self.frStatus)
        self.btnExport.grid(column=2, row=0, sticky='e', padx=5, pady=6)
        self.apply_style(self.btnExport)
        self.btnExport.configure(text='''Export table''')
        self.btnExport.configure(width=20)
        self.btnExport.configure(height=2)
        self.btnExport.configure(command=lambda: self.export_dialog())

        self.cbShowCalc = Checkbutton(self.frStatus)
        self.cbShowCalc.grid(row=0, column=5, pady=5, padx = 5, sticky='ew')
        self.apply_style(self.cbShowCalc)
        self.cbShowCalc.configure(justify=LEFT)
        self.cbShowCalc.configure(text='''Show calc?''')
        self.cbShowCalc.configure(variable=gui_support.varShowCalc)
        self.cbShowCalc.configure(command=lambda: self.plot_text(self.set_pval_dval()[0], self.set_pval_dval()[1]))

        # _____________________frGraph___________________________________________________________________________________
        self.frGraph = Frame(master)
        self.frGraph.configure(relief=GROOVE)
        self.frGraph.configure(borderwidth="2")
        self.frGraph.configure(relief=GROOVE)
        self.frGraph.configure(background="#d9d9d9")
        self.frGraph.configure(highlightbackground="#d9d9d9")
        self.frGraph.configure(highlightcolor="black")
        self.frGraph.grid(row=0, column=1, rowspan=3, sticky='nswe')
        self.frGraph.columnconfigure(0, weight=1)
        self.frGraph.rowconfigure(0, weight=1)
        self.frGraph.rowconfigure(1, weight=0)
        self.frGraph.rowconfigure(2, weight=1)
        self.frGraph.rowconfigure(3, weight=0)
        self.frGraph.rowconfigure(4, weight=1)
        self.frGraph.rowconfigure(5, weight=0)

        #______________frCon
        self.frConc = Frame(self.frGraph)
        self.frConc.grid(row=0, sticky='nswe')
        self.frConc.configure(relief=GROOVE)
        self.frConc.configure(borderwidth="2")
        self.frConc.configure(relief=GROOVE)
        self.frConc.configure(background="#d9d9d9")
        self.frConc.configure(highlightbackground="#d9d9d9")
        self.frConc.configure(highlightcolor="black")

        self.fig = Figure(figsize=(6, 2.15), frameon=False)
        self.ax_conc = self.fig.add_subplot(111)
        self.ax_conc.axes
        self.ax_conc.set_xlabel('207Pb/235U')
        self.ax_conc.set_ylabel('206Pb/238U')
        self.ax_conc.set_title('Concordia')
        try:
            self.ax_conc.plot(list(range(0, EarthAge)), graph_to_draw)
        except UnboundLocalError:
            pass

        self.canvas_conc = FigureCanvasTkAgg(self.fig, self.frConc)
        self.canvas_conc.draw()
        self.canvas_conc.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.frConcToolbar = Frame(self.frGraph)
        self.frConcToolbar.grid(row=1, sticky='ew')
        self.frConcToolbar.configure(relief=GROOVE)
        self.frConcToolbar.configure(borderwidth="2")
        self.frConcToolbar.configure(relief=GROOVE)
        self.frConcToolbar.configure(background="#d9d9d9")
        self.frConcToolbar.configure(highlightbackground="#d9d9d9")
        self.frConcToolbar.configure(highlightcolor="black")
        self.frConcToolbar.configure(width=200)

        # ______________frProb
        self.frProb = Frame(self.frGraph)
        self.frProb.grid(row=2, sticky='nswe')
        self.frProb.configure(relief=GROOVE)
        self.frProb.configure(borderwidth="2")
        self.frProb.configure(relief=GROOVE)
        self.frProb.configure(background="#d9d9d9")
        self.frProb.configure(highlightbackground="#d9d9d9")
        self.frProb.configure(highlightcolor="black")

        try:
            if g_graph_settings.pdp_kde_hist == 0:
                graph_to_draw = g_grainset.kde(g_graph_settings.bandwidth)

            elif g_graph_settings.pdp_kde_hist == 1:
                graph_to_draw = g_grainset.pdp()[0]
        except NameError:
            pass

        self.fig = Figure(figsize=(6, 2.15), frameon=False)
        self.ax_prob = self.fig.add_subplot(111)
        self.ax_prob.axes
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
        self.frProbToolbar.grid(row=3, sticky='ew')
        self.frProbToolbar.configure(relief=GROOVE)
        self.frProbToolbar.configure(borderwidth="2")
        self.frProbToolbar.configure(relief=GROOVE)
        self.frProbToolbar.configure(background="#d9d9d9")
        self.frProbToolbar.configure(highlightbackground="#d9d9d9")
        self.frProbToolbar.configure(highlightcolor="black")
        self.frProbToolbar.configure(width=200)

        # ______________frCum
        self.frCum = Frame(self.frGraph)
        self.frCum.grid(row=4, sticky='nswe')
        self.frCum.configure(relief=GROOVE)
        self.frCum.configure(borderwidth="2")
        self.frCum.configure(relief=GROOVE)
        self.frCum.configure(background="#d9d9d9")
        self.frCum.configure(highlightbackground="#d9d9d9")
        self.frCum.configure(highlightcolor="black")
        self.frCum.configure(height=10)

        try:
            if g_graph_settings.pdp_kde_hist == 0:
                graph_to_draw = g_grainset.ckde(g_graph_settings.bandwidth)

            elif g_graph_settings.pdp_kde_hist == 1:
                graph_to_draw = g_grainset.cpdp()
        except NameError:
            pass

        self.fig = Figure(figsize=(6, 2.15), frameon=False)
        self.ax_cum = self.fig.add_subplot(111)
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
        self.frCumToolbar.grid(row=5, sticky='ew')
        self.frCumToolbar.configure(relief=GROOVE)
        self.frCumToolbar.configure(borderwidth="2")
        self.frCumToolbar.configure(relief=GROOVE)
        self.frCumToolbar.configure(background="#d9d9d9")
        self.frCumToolbar.configure(highlightbackground="#d9d9d9")
        self.frCumToolbar.configure(highlightcolor="black")
        self.frCumToolbar.configure(width=200)

        global toolbarConc, toolbarProb, toolbarCum
        toolbarConc = NavigationToolbar2Tk(self.canvas_conc, self.frConcToolbar)
        toolbarProb = NavigationToolbar2Tk(self.canvas_prob, self.frProbToolbar)
        toolbarCum = NavigationToolbar2Tk(self.canvas_cum, self.frCumToolbar)

        #________________Menu___________________________________________________________________________________________
        self.menubar = Menu(master, font="TkMenuFont", bg=_bgcolor, fg=_fgcolor)
        master.configure(menu=self.menubar)
        mFile = Menu(self.menubar, tearoff=False)
        mEdit = Menu(self.menubar, tearoff=False)
        mAbout = Menu(self.menubar, tearoff=False)

        self.menubar.add_cascade(label="File", menu=mFile, underline=0)
        self.menubar.add_cascade(label="Edit", menu=mEdit, underline=0)
        self.menubar.add_cascade(label="About", menu=mAbout, underline=0)

        mFile.add_command(label="New Session", underline=0, accelerator="Ctrl+N",
                          command=lambda: self.reset_controls(False))


        mFile.add_command(label="Open Session", underline=0, accelerator="Ctrl+O",
                          command=lambda: filedialog.askopenfilename(initialdir="/", title="Select File",
                                                                     filetypes=(("Text files", "*.txt"),
                                                              ("Comma separated values files", "*.csv"),
                                                                                 ("All files", "*.*"))))

        mFile.add_command(label="Save Session", underline=0, accelerator="Ctrl+S", command=lambda: print("Save"))

        mFile.add_separator()

        mFile.add_command(label="Import Data", underline=0, accelerator="Ctrl+I", command=lambda:
                                  self.open_and_load_file(g_filters, self.Table, g_list_col_names))


        mFile.add_command(label="Export Table", underline=7, accelerator="Ctrl+E+T")
        mFile.add_command(label="Export Graph", underline=7, accelerator="Ctrl+E+G")
        mFile.add_separator()
        mFile.add_command(label="Exit", underline=1, command=root.quit, accelerator="Ctrl+Q")
        mEdit.add_command(label="Undo", accelerator="Ctrl+Z")
        mEdit.add_command(label="Redo", accelerator="Ctrl+Shift+Z")
        mEdit.add_separator()
        mEdit.add_command(label="Settings")

        self.reset_controls(False)

    #_____________Class Methods_________________________________________________________________________________________
    def apply_style(self, obj):
        obj.configure(activebackground="#f9f9f9")
        obj.configure(activeforeground="black")
        obj.configure(background="#d9d9d9")
        obj.configure(disabledforeground="#a3a3a3")
        obj.configure(foreground="#000000")
        obj.configure(highlightbackground="#d9d9d9")
        obj.configure(highlightcolor="black")

    def open_and_load_file(self):
        try:
            try:
                global g_plot_txt, g_directory, g_file_type, g_filters, g_list_col_names
                if g_plot_txt != "":
                    g_plot_txt.remove()
                g_filters.sample_name_filter = []
                user_file = filedialog.askopenfilename(
                    initialdir=g_directory, title="Select file", filetypes=(("Text files", "*.txt"),
                                                                    ("Comma separated values files", "*.csv"),
                                                                    ("All files", "*.*")))

                g_directory = os.path.split(user_file)[0]
                root.title(user_file + ' — De-Zir-teer')
                an_set = []
                file = imported_file(user_file)
                g_file_type = file[1]
                for i in range(1, file[2]):
                    an = file_to_analysis(file, i)
                    if str(an.analysis_name) == "TS16_20":
                        pass
                    an_set.append(an)

                global g_grainset
                g_grainset = AnalysesSet(an_set, 'set#1')
                g_grainset.good_bad_sets(g_filters)

                global pars_onChange
                pars_onChange = [g_filters, self.Table, g_grainset, g_list_col_names]

                sys.stdout.flush()
                global g_number_of_good_grains
                g_number_of_good_grains = gui_support.fill_data_table(self.Table, g_grainset, g_filters, g_list_col_names)
                global g_list_of_samples
                g_list_of_samples = same_sample_set(g_grainset, str(self.cbSeparatorType.get()))
                self.reset_controls(True)
                self.clear_prev_or_remove_text()

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
        global g_plot_txt, g_prev_n, g_prev_cum
        g_prev_n = 0
        g_prev_cum = []
        self.ax_conc.clear()
        self.ax_prob.clear()
        self.ax_cum.clear()
        self.canvas_conc.draw()
        self.canvas_cum.draw()
        self.canvas_prob.draw()
        self.btnClear.configure(state=DISABLED)
        g_plot_txt = ""

    def export_dialog(self):
        file_main = filedialog.asksaveasfile(mode='w', defaultextension=".csv", initialdir=g_directory,
                                             filetypes=(("Comma separated values files", "*.csv"),
                                                          ("All files", "*.*")))
        file_prob = os.path.dirname(str(file_main.name)) + '/' + \
                    os.path.splitext(os.path.basename(str(file_main.name)))[0]+'_prob_cum' + '.csv'

        gui_support.export_table(g_grainset, g_filters, g_list_col_names, g_graph_settings, file_main, file_prob)

    def reset_controls(self, is_data_present):
        unfinished_features = [self.chbAnchored, self.entAnchoredAge, self.chbFitDiscordia, self.rbUseCorr,
                               self.cbTypePbc, self.rbDiscUbased, self.rbUBased,
                               self.chbInclude207235Err, self.scErrFilter, self.scUconcCutoff]
        if is_data_present:
            for var_frame in (self.frImport, self.frFilter, self.frDisc, self.frGraphSettings, self.frStatus):
                for child in var_frame.winfo_children():
                    if child not in unfinished_features:
                        child.configure(state=NORMAL)
            self.rbDiscUbased.configure(state=DISABLED)
            self.rbDiscAgeFixedLim.configure(state=DISABLED)
            self.rbDisc67_68.configure(state=DISABLED)
            self.rbDisc75_68.configure(state=DISABLED)
            self.cbSeparatorType.configure(state="readonly")
            self.rbAgeFixedLim.select()
            self.chbDiscLinked2Age.select()
            self.rbDiscAgeFixedLim.select()

            self.lboxSamples.delete(0, END)
            for item in g_list_of_samples:
                self.lboxSamples.insert(END, item.name)
            if self.lboxSamples.get(0) == '':
                status_text = ' data, bad divider'
                status_color = 'red'
            else:
                status_text = ' data OK'
                status_color = 'green'
            self.lbShowStatus.configure(text=g_file_type+status_text, fg=status_color)
        else:
            self.lboxSamples.delete(0, END)
            for var_frame in (self.frImport, self.frFilter, self.frDisc, self.frGraphSettings, self.frStatus):
                for child in var_frame.winfo_children():
                    child.configure(state=DISABLED)
            self.btnImport.configure(state='normal')
            self.lbSeparatorType.configure(state='normal')
            self.cbSeparatorType.configure(state='readonly')
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

    #adds or removes text to the cum_plot, depending on the checked state of the cbShowCalc
    def plot_text(self, pval, dval):
        global g_plot_txt
        if g_graph_settings.pdp_kde_hist == 0:
            peaks = g_grainset.kde(g_graph_settings.bandwidth)[1]
        else:
            peaks=g_grainset.pdp()[1]
        if gui_support.varShowCalc.get() == 1:
            text_to_show = \
                "n={}\n" \
                "Min age={}; " \
                "Max age={}\n" \
                "WA age={}±{}(2σ int.);" \
                "±{}(95%conf)\n" \
                "MSWD={}\n" \
                "KS p-value={}; " \
                "d-value={}\n" \
                "peaks at={}".format(
                    g_number_of_good_grains[0],
                    int(g_number_of_good_grains[6]),
                    int(g_number_of_good_grains[5]),
                    round((g_number_of_good_grains[1]), 1),
                    2 * round((g_number_of_good_grains[2]), 1),
                    round((g_number_of_good_grains[3]), 1),
                    int(g_number_of_good_grains[4]),
                    round(pval, 2),
                    round(dval, 2),
                    peaks
                )


        else:
            if g_plot_txt != "":
                g_plot_txt.remove()
            text_to_show = ""

        g_plot_txt = self.ax_cum.text(0.05, 0.40, text_to_show, transform=self.ax_cum.transAxes)
        self.plot_peaks()
        self.canvas_cum.draw()
        self.canvas_prob.draw()

    #calculates KS p- and d-values from current and previous grainsets. If previous grainset doesn't exist,
    #sets 0 values to both variables
    def set_pval_dval(self):
        global g_prev_cum, g_grainset
        if g_prev_cum == []:
            pval = 0
            dval = 0
        else:
            if g_graph_settings.pdp_kde_hist == 0:
                curr_cum = g_grainset.ckde(g_graph_settings.bandwidth)
            elif g_graph_settings.pdp_kde_hist == 1:
                curr_cum = g_grainset.cpdp()
            else:
                curr_cum = []
            dval = d_value(curr_cum, g_prev_cum)
            pval = p_value(dval, g_number_of_good_grains[0], g_prev_n[0])
        return [pval, dval]


    def min_max_ages(self):
        # choosing age interval based on user's input
        if gui_support.varLimitAgeSpectrum.get() == 1:
            min_age = g_number_of_good_grains[6]
            max_age = g_number_of_good_grains[5]
        else:
            min_age = 1
            max_age = EarthAge
        return [min_age, max_age]

    def concordia_type(self):
        # choosing concordia type base on user's input
        if g_graph_settings.conc_type == 0:  # if conventional concordia
            conc_graph_x = [i[1] for i in concordia_table]
            conc_graph_y = [i[0] for i in concordia_table]
            conc_title = "Conventional Concordia"
            conc_graph_xtitle = "207Pb/235U"
            conc_graph_ytitle = "206Pb/238U"
            xconc = 1
            yconc = 0
        else:  # Tera-Wasserburgh
            conc_graph_x = [(1 / i[0]) for i in concordia_table]
            conc_graph_y = [i[2] for i in concordia_table]
            conc_title = "Tera-Wasserburg Concordia"
            conc_graph_xtitle = "238U/206Pb"
            conc_graph_ytitle = "207Pb/206Pb"
            xconc = 3
            yconc = 2
        return [conc_graph_x, conc_graph_y, conc_title, conc_graph_xtitle, conc_graph_ytitle, xconc, yconc]

    def kde_pdp_hist(self):
        # choosing kde/pdp/hist based on user input
        if g_graph_settings.pdp_kde_hist == 0:
            prob_graph_to_draw = g_grainset.kde(g_graph_settings.bandwidth)[0]
            cum_graph_to_draw = g_grainset.ckde(g_graph_settings.bandwidth)
            prob_title = "Kernel Density Estimates (KDE)"
            cum_title = "Cumulative KDE"
        elif g_graph_settings.pdp_kde_hist == 1:
            prob_graph_to_draw = g_grainset.pdp()[0]
            cum_graph_to_draw = g_grainset.cpdp()
            prob_title = "Probability Density Plot (PDP)"
            cum_title = "Cumulative PDP"
        else:
            tuple_list = sorted(list(g_grainset.good_set.values()), key=lambda x: x[0])
            prob_graph_to_draw = [x[0] for x in tuple_list]
            cum_graph_to_draw=[]
            prob_title = "Histogram"
            cum_title = "Cumulative Histogram"
        return[prob_graph_to_draw, cum_graph_to_draw, prob_title, cum_title]


    def draw_concordia_ticks(self, xconc, yconc):
        for t in [100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500]:
            x = calc_ratio(t)[xconc]
            y = calc_ratio(t)[yconc]
            self.ax_conc.plot(x,y,'ko')
            self.ax_conc.text(x, y, str(t), style='italic', )
            # self.ax_conc.annotate(str(t), xy=(x, y), xytext=(x-0.4*x, y+0.4*y),
            #            arrowprops=dict(facecolor='black', shrink=0.05))

    def plot_conc_ellipses(self, args):
        # plots ellipses on concordia-discordia diagram
        for zir in g_grainset.good_set:
            sigma_level = g_graph_settings.eclipses_at

            # conventional concordia
            if g_graph_settings.conc_type == 0:
                corr_coef = zir.corr_coef_75_68
                x_conc = zir.pb207_u235[0]  # x-center of the oval
                y_conc = zir.pb206_u238[0]  # y-center of the oval
                x_err = zir.pb207_u235[1]
                y_err = zir.pb206_u238[1]
            # Tera-Wasserburg concordia
            else:
                corr_coef = zir.corr_coef_86_76
                j = zir.u238_pb206()
                x_conc = j[0]
                x_err = j[1]
                y_conc = zir.pb207_pb206[0]
                y_err = zir.pb207_pb206[1]

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

            if args != "":
                if zir.analysis_name == args[0]:
                    oval_color = 'blue'
                    oval_fill = True
                else:
                    oval_color = 'red'
                    oval_fill = False
            else:
                oval_color = 'red'
                oval_fill = False
            el = Ellipse(xy=(x_conc, y_conc), width=a * 2, height=b * 2, angle=degrees(ang), color=oval_color,
                         fill=oval_fill)
            self.ax_conc.add_patch(el)

    def plot_hist(self, min_age, max_age, prob_graph_to_draw):
        bin_sequence = []
        age = min_age
        bin_width = self.scBinwidth.get()
        while age < max_age:
            bin_sequence.append(age)
            age += bin_width
        self.ax_prob.hist(prob_graph_to_draw, bins=bin_sequence, density=True, cumulative=False)
        self.ax_cum.hist(prob_graph_to_draw, bins=bin_sequence, density=True, cumulative=True)

    def set_axes(self, conc_title,conc_graph_xtitle, conc_graph_ytitle, prob_title, cum_title, conc_graph_x,
                 conc_graph_y, min_age, max_age):
        # set axis of all graphs
        self.ax_conc.set_title(conc_title)
        self.ax_conc.set_xlabel(conc_graph_xtitle, labelpad=-12, fontsize=8, position=(0.54, 1e6))
        self.ax_conc.set_ylabel(conc_graph_ytitle, fontsize=8)
        self.ax_prob.set_title(prob_title)
        self.ax_prob.set_xlabel('Age (Ma)', labelpad=-12, fontsize=8, position=(0.54, 1e6))
        self.ax_cum.set_title(cum_title)
        self.ax_cum.set_xlabel('Age (Ma)', labelpad=-12, fontsize=8, position=(0.54, 1e6))
        self.ax_conc.plot(conc_graph_x[min_age: max_age], conc_graph_y[min_age: max_age])

    def plot_peaks(self):
        prob_graph_to_draw = self.kde_pdp_hist()[0]
        min_max_age = self.min_max_ages()
        min_age = min_max_age[0]
        max_age = min_max_age[1]
        self.ax_prob.clear()
        self.canvas_prob.draw()
        self.ax_prob.plot(list(range(min_age, max_age)), prob_graph_to_draw[min_age: max_age])
        if gui_support.varShowCalc.get() == 1:
            i = 0
            if g_graph_settings.pdp_kde_hist == 0:
                list_peaks = g_grainset.kde(g_graph_settings.bandwidth)[1]
            elif g_graph_settings.pdp_kde_hist == 1:
                list_peaks = g_grainset.pdp()[1]
            else:
                list_peaks = []
            while i<len(list_peaks):
                self.ax_prob.axvline(list_peaks[i], color='black')
                i += 1
        else:
           pass

    def prob_cum_plot(self, min_age, max_age, prob_graph_to_draw, cum_graph_to_draw):
        self.ax_cum.plot(list(range(min_age, max_age)), cum_graph_to_draw[min_age: max_age])
        self.plot_peaks() #ax_prob.plot is done here


    def prob_cum_hist_plot(self, do_hist, min_age, max_age, prob_graph_to_draw, cum_graph_to_draw):
        if not do_hist:
            self.prob_cum_plot(min_age, max_age, prob_graph_to_draw, cum_graph_to_draw)
        else:
            self.plot_hist(min_age, max_age, prob_graph_to_draw)

    def clear_prev_or_remove_text(self):
        # clears previous graph, if user chooses to in the cbKeepPrev, else just removes text from cum_plot
        global g_plot_txt
        if gui_support.varKeepPrev.get() == 0:
            self.clear_graph()
        else:
            if g_plot_txt != "":
                g_plot_txt.remove()
        g_plot_txt = ""

    def plot_conc_text_peaks(self):
        global g_prev_n, g_prev_cum
        self.plot_text(self.set_pval_dval()[0], self.set_pval_dval()[1])
        self.canvas_conc.draw()
        #self.canvas_prob.draw() and self.canvas_cum.draw are executed in plot_text
        self.btnClear.configure(state=NORMAL)
        g_prev_n = g_number_of_good_grains
        if g_graph_settings.pdp_kde_hist == 0:
            g_prev_cum = g_grainset.ckde(g_graph_settings.bandwidth)
        else:
            g_prev_cum = g_grainset.cpdp()

    #draws the graph based on the data and user settings. Clears the previous graph, or draws on top of it,
    #depending on user settings
    def clear_and_plot(self, *args):
        global g_filters, g_grainset, g_number_of_good_grains, g_plot_txt, g_prev_cum, g_prev_n
        g_filters.sample_name_filter = []

        #gets the user-selected items from the listbox
        item_indexes = self.lboxSamples.curselection()
        items = [self.lboxSamples.get(item_indexes) for item_indexes in item_indexes]
        g_filters.sample_name_filter = items
        g_number_of_good_grains = gui_support.fill_data_table(self.Table, g_grainset, g_filters, g_list_col_names)

        #checks if histogram is to be drawn
        do_hist = (g_graph_settings.pdp_kde_hist == 2)

        # cropping age interval: either full, or cropped from min_age to max_age
        age_lim = self.min_max_ages()
        min_age = age_lim[0]
        max_age = age_lim[1]

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
        prob_graph_to_draw = l_kde_pdp_hist[0]
        cum_graph_to_draw = l_kde_pdp_hist[1]
        prob_title = l_kde_pdp_hist[2]
        cum_title = l_kde_pdp_hist[3]

        # set axis of all graphs
        self.set_axes(conc_title, conc_graph_xtitle, conc_graph_ytitle, prob_title, cum_title, conc_graph_x,
                      conc_graph_y, min_age, max_age)

        self.draw_concordia_ticks(xconc, yconc)

        if args:
            user_selected_analysis = args
        else:
            user_selected_analysis = ""
        try:
            #plots ellipses on concordia-discordia diagram
            self.plot_conc_ellipses(user_selected_analysis)

            # plotting KDE/CKDE, PDP/CPDP or histogram
            self.prob_cum_hist_plot(do_hist, min_age, max_age, prob_graph_to_draw, cum_graph_to_draw)

        except ValueError:
            pass

        except TypeError:
            pass

        finally:
            self.plot_conc_text_peaks()



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


def main():
    global root, g_list_col_names, g_grainset, g_filters, g_graph_settings, prob_fig, prob_subplot
    global g_list_of_samples, g_directory, g_number_of_good_grains, g_prev_cum, g_prev_n
    g_prev_cum = []
    g_directory = "C:\odrive\Amazon Cloud Drive\cloud\Geochron\Santa Cruz LA"
    g_list_col_names = ['208Pb/232Th', '208/232±1s',
                        '207Pb/206Pb', '207/206±1s',
                        '207Pb/235U', '207/235±1s',
                        '206Pb/238U', '206/238±1s',
                        'corr. coef.75_68', 'corr. coef.86_76',
                        'Uconc (approx. ppm)', 'Uconc±1s',
                        'pbc (approx. ppm)', 'pbc±1s',
                        '206Pb/204Pb', '206Pb/204Pb±1s',
                        '207Pb/204Pb', '207Pb/204Pb±1s',
                        '208Pb/204Pb', '208Pb/204Pb±1s',
                        'Age 208U/232Th', 'Age208/232±1s',
                        'Age 207Pb/206Pb', 'Age207/206±1s',
                        'Age 207Pb/235U', 'Age207/235±1s',
                        'Age 206Pb/238U', 'Age206/238±1s',
                        'disc. 207/206-206/238', 'disc. 207/235-206/238',
                        'is grain good?', 'best age system', 'best age', 'best age±1s']
    fill_pbpb_table()
    fill_concordia_table()
    g_filters = Filters()
    g_prev_n = 0
    g_graph_settings = gui_support.GraphSettings()
    root = Tk()
    root.title('De-Zir-teer')
    root.wm_resizable(1, 1)
    gui_support.set_Tk_var()
    master = OperationWindow(root)
    if __name__ == "__main__":
        root.mainloop()


main()
