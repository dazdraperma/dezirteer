import unittest
import sys
sys.path.append('../dezirteer')
import math_module
import main_gui
import gui_support

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


class TestMath(unittest.TestCase):

    def set_vars(self):
        global g_filters, g_table, g_grainset, g_list_col_names, g_operwindow, g_pars_onChange
        g_operwindow = main_gui.OperationWindow(Tk())
        g_filters = main_gui.g_filters
        g_table = g_operwindow.Table
        g_grainset = main_gui.g_grainset
        g_list_col_names = main_gui.g_list_col_names
        g_pars_onChange = [g_filters, g_table, g_grainset, g_list_col_names]


    def setUp(self):
        global g_filters, g_table, g_grainset, g_list_col_names, g_operwindow, g_pars_onChange
        self.zir = math_module.Analysis('test_zircon', 15, (0.2003, 0.0008, 0.0046), (2.082, 0.009, 0.07), 0.6, 0.6,
                                   (0.0617, 0.0003, 0.0003), (0.758, 0.0003, 0.0015), (0, 0, 0), (0, 0, 0), (0, 0, 0),
                                   (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), 1)

        self.set_vars()


    def test_calc_ratio(self):
        result = math_module.calc_ratio(1000)[0]
        self.assertEqual(result, 0.16780392747297124)

    def test_analysis(self):
        self.assertEqual(round(self.zir.calc_age(0)[0], 0), 1177)

    def test_import_plot(self):
        global g_filters, g_table, g_grainset, g_list_col_names, g_operwindow, g_pars_onChange

        g_operwindow.open_and_load_file("C:/Program Files (x86)/Dezirteer/Examples/iolite_example.txt")
        g_operwindow.entPosDiscFilt.delete(0, "end")
        g_operwindow.entPosDiscFilt.insert(0, 4)


        self.set_vars()

        gui_support.onChange(6, 4, g_pars_onChange)
        gui_support.onChange(7, 10, g_pars_onChange)
        g_operwindow.clear_and_plot()
        self.assertEqual(main_gui.g_number_of_good_grains[0], 165)
        self.assertEqual(int(main_gui.g_number_of_good_grains[1]), 877)
        self.assertEqual(round(main_gui.g_number_of_good_grains[2], 2), 0.39)
        self.assertEqual(int(main_gui.g_number_of_good_grains[3]), 114)
        self.assertEqual(int(main_gui.g_number_of_good_grains[4]), 22555)
        self.assertEqual(main_gui.g_number_of_good_grains[5], 3265)
        self.assertEqual(main_gui.g_number_of_good_grains[6], 404)


        gui_support.onChange(6, 20, g_pars_onChange)
        g_operwindow.clear_and_plot()
        self.assertEqual(main_gui.g_number_of_good_grains[0], 173)
        self.assertEqual(int(main_gui.g_number_of_good_grains[1]), 908)
        self.assertEqual(round(main_gui.g_number_of_good_grains[2], 2), 0.38)
        self.assertEqual(int(main_gui.g_number_of_good_grains[3]), 116)
        self.assertEqual(int(main_gui.g_number_of_good_grains[4]), 23593)
        self.assertEqual(main_gui.g_number_of_good_grains[5], 3265)
        self.assertEqual(main_gui.g_number_of_good_grains[6], 404)



if __name__ == '__main__':
    unittest.main()