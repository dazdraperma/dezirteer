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


class TestGui(unittest.TestCase):

    def set_vars(self):
        global g_filters, g_table, g_grainset, g_list_col_names, g_operwindow, g_pars_onChange, g_file_location_and_name, \
            g_matrix_results, g_matrix_inputs

        g_operwindow = main_gui.OperationWindow(Tk())
        g_filters = main_gui.g_filters
        g_table = g_operwindow.Table
        g_grainset = main_gui.g_grainset
        g_list_col_names = main_gui.g_list_col_names
        g_pars_onChange = [g_filters, g_table, g_grainset, g_list_col_names]
        g_file_location_and_name = "iolite_example.txt"
        g_matrix_inputs = [[4, 10], [20, 10]]
        g_matrix_results = [[165, 877, 0.39, 114, 22555, 3265, 404], [173, 908, 0.38, 116, 23593, 3265, 404]]

    '''def verify_values(self, n, wa, sigma, ninety_five, mswd, max_age, min_age):
        is_n_good = (main_gui.g_number_of_good_grains[0] ==  n)
        is_wa_good = (int(main_gui.g_number_of_good_grains[1]) == wa)
        is_sigma_good = (round(main_gui.g_number_of_good_grains[2], 2) == sigma)
        is_ninety_five_good = (int(main_gui.g_number_of_good_grains[3]) == ninety_five)
        is_mswd_good = (int(main_gui.g_number_of_good_grains[4]) == mswd)
        is_max_age_good = (main_gui.g_number_of_good_grains[5] == max_age)
        is_min_age_good = (main_gui.g_number_of_good_grains[6] == min_age)
        return [is_n_good, is_wa_good, is_sigma_good, is_ninety_five_good, is_mswd_good, is_max_age_good, is_min_age_good]'''

    def setUp(self):
        self.zir = math_module.Analysis('test_zircon', 15, (0.2003, 0.0008, 0.0046), (2.082, 0.009, 0.07), 0.6, 0.6,
                                   (0.0617, 0.0003, 0.0003), (0.758, 0.0003, 0.0015), (0, 0, 0), (0, 0, 0), (0, 0, 0),
                                   (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), 1)
        self.set_vars()

    def test_calc_ratio(self):
        result = math_module.calc_ratio(1000)[0]
        self.assertEqual(result, 0.16780392747297124)

    def test_analysis(self):
        self.assertEqual(round(self.zir.calc_age(0)[0], 0), 1177)

    def test_discordances(self):
        g_operwindow.open_and_load_file(g_file_location_and_name)
        self.set_vars()
        gui_support.onChange(6, g_matrix_inputs[0][0], g_pars_onChange)
        gui_support.onChange(7, g_matrix_inputs[0][1], g_pars_onChange)
        g_operwindow.clear_and_plot()

        #ver_values = self.verify_values(165, 877, 0.39, 114, 22555, 3265, 404)
        #this is the shortest code, but the most difficult to debug
        '''for i in ver_values:
            self.assertEqual(ver_values[i], True)
            print (i)'''
        # this is a longer code, but less difficult to debug
        '''self.assertEqual(ver_values[0], True)
        self.assertEqual(ver_values[1], True)
        self.assertEqual(ver_values[2], True)
        self.assertEqual(ver_values[3], True)
        self.assertEqual(ver_values[4], True)
        self.assertEqual(ver_values[5], True)
        self.assertEqual(ver_values[6], True)'''

        # this is the longest and least neat code, but very easy to debug
        self.assertEqual(main_gui.g_number_of_good_grains[0], g_matrix_results[0][0])
        self.assertEqual(int(main_gui.g_number_of_good_grains[1]), g_matrix_results[0][1])
        self.assertEqual(round(main_gui.g_number_of_good_grains[2], 2), g_matrix_results[0][2])
        self.assertEqual(int(main_gui.g_number_of_good_grains[3]), g_matrix_results[0][3])
        self.assertEqual(int(main_gui.g_number_of_good_grains[4]), g_matrix_results[0][4])
        self.assertEqual(main_gui.g_number_of_good_grains[5], g_matrix_results[0][5])
        self.assertEqual(main_gui.g_number_of_good_grains[6], g_matrix_results[0][6])


        gui_support.onChange(6, g_matrix_inputs[1][0], g_pars_onChange)
        gui_support.onChange(7, g_matrix_inputs[1][1], g_pars_onChange)

        g_operwindow.clear_and_plot()
        self.assertEqual(main_gui.g_number_of_good_grains[0], g_matrix_results[1][0])
        self.assertEqual(int(main_gui.g_number_of_good_grains[1]), g_matrix_results[1][1])
        self.assertEqual(round(main_gui.g_number_of_good_grains[2], 2), g_matrix_results[1][2])
        self.assertEqual(int(main_gui.g_number_of_good_grains[3]), g_matrix_results[1][3])
        self.assertEqual(int(main_gui.g_number_of_good_grains[4]), g_matrix_results[1][4])
        self.assertEqual(main_gui.g_number_of_good_grains[5], g_matrix_results[1][5])
        self.assertEqual(main_gui.g_number_of_good_grains[6], g_matrix_results[1][6])

        gui_support.onChange(6, 1, g_pars_onChange)
        gui_support.onChange(7, 1, g_pars_onChange)
        g_operwindow.clear_and_plot()
        self.assertEqual(main_gui.g_number_of_good_grains[0], 67)
        self.assertEqual(int(main_gui.g_number_of_good_grains[1]), 830)
        self.assertEqual(round(main_gui.g_number_of_good_grains[2], 2), 0.64)
        self.assertEqual(int(main_gui.g_number_of_good_grains[3]), 171)
        self.assertEqual(int(main_gui.g_number_of_good_grains[4]), 18181)
        self.assertEqual(main_gui.g_number_of_good_grains[5], 2930)
        self.assertEqual(main_gui.g_number_of_good_grains[6], 408)


if __name__ == '__main__':
    unittest.main()