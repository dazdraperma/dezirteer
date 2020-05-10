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
    global g_operwindow
    g_operwindow = main_gui.OperationWindow(Tk())

    def set_vars(self):
        global g_filters, g_table, g_grainset, g_list_col_names, g_operwindow, g_pars_onChange, \
            g_matrix_filenames, g_matrix_results, g_matrix_inputs

        #g_operwindow = main_gui.OperationWindow(Tk())
        g_filters = main_gui.g_filters
        g_table = g_operwindow.Table
        g_grainset = main_gui.g_grainset
        g_list_col_names = main_gui.g_list_col_names
        g_pars_onChange = [g_filters, g_table, g_grainset, g_list_col_names]

        g_matrix_filenames = ["iolite_example.txt",
                              "glitter_example.txt"]

        #Format:
        #positive_disc,
        #negative_disc (omit minus sign),
        #use err filter (0==no, 1==yes),
        #err filter value,
        #207/235 used (0==no, 1==yes)
        g_matrix_inputs = [
            [4, 10, 0, 5, 1],   #test 0
            [20, 10, 0, 5, 0],  #test 1
            [1, 1, 0, 5, 0],    #test 2
            [20, 10, 1, 1.2, 0],#test 3
            [20, 10, 1, 2.5, 1] #test 4
        ]

        # Format:
        # file0[n, WA, 1sigma, 95%, MSWD, maxage, minage],
        # File1[n, WA, 1sigma, 95%, MSWD, maxage, minage],
        # ...
        # FileN[n, WA, 1sigma, 95%, MSWD, maxage, minage]
        g_matrix_results = [
            [#file 0
                [165, 877, 0.39, 114, 22555, 3265, 404],#file 0, test 0
                [173, 908, 0.38, 116, 23593, 3265, 404],#file 0, test 1
                [67, 830, 0.64, 171, 18181, 2930, 408], #file 0, test 2
                [169, 908, 0.38, 117, 24145, 3265, 404],#file 0, test 3
                [122, 725, 0.41, 97, 14388, 2665, 404]  #file 0, test 4
            ],
            [#file 1
                [104, 267, 0.4, 31, 1604, 1081, 133], #file 1, test 0
                [146, 258, 0.32, 22, 1289, 2695, 133],#file 1, test 1
                [44, 311, 0.73, 63, 1895, 1081, 177], #file 1, test 2
                [13, 236, 0.85, 69, 1406, 2695, 160], #file 1, test 3
                [15, 260, 0.87, 62, 1113, 767, 160]   #file 1, test 4
            ]
        ]

    def setUp(self):
        self.set_vars()

    def test_matrices_file(self):
        for file_num in range(len(g_matrix_filenames)):
            g_operwindow.open_and_load_file(g_matrix_filenames[file_num], False)
            self.set_vars()

            for test_num in range(len(g_matrix_inputs)):
                gui_support.onChange(6, g_matrix_inputs[test_num][0], g_pars_onChange)
                gui_support.onChange(7, g_matrix_inputs[test_num][1], g_pars_onChange)

                gui_support.onChange(5, g_matrix_inputs[test_num][2], g_pars_onChange, g_matrix_inputs[test_num][3],
                                     g_operwindow.entErrFilter, g_operwindow.chbInclude207235Err)
                gui_support.onChange(20,  g_matrix_inputs[test_num][3], g_pars_onChange)
                gui_support.onChange(22, g_matrix_inputs[test_num][4], g_pars_onChange)

                g_operwindow.clear_and_plot()
                self.assertEqual(main_gui.g_number_of_good_grains[0], g_matrix_results[file_num][test_num][0])
                self.assertEqual(int(main_gui.g_number_of_good_grains[1]), g_matrix_results[file_num][test_num][1])
                self.assertEqual(round(main_gui.g_number_of_good_grains[2], 2), g_matrix_results[file_num][test_num][2])
                self.assertEqual(int(main_gui.g_number_of_good_grains[3]), g_matrix_results[file_num][test_num][3])
                self.assertEqual(int(main_gui.g_number_of_good_grains[4]), g_matrix_results[file_num][test_num][4])
                self.assertEqual(main_gui.g_number_of_good_grains[5], g_matrix_results[file_num][test_num][5])
                self.assertEqual(main_gui.g_number_of_good_grains[6], g_matrix_results[file_num][test_num][6])


if __name__ == '__main__':
    unittest.main()