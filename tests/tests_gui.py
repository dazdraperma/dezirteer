import unittest
import sys
sys.path.append('../dezirteer')
import math_module
import main_gui
import gui_support
import datetime

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
        #how to calc best age (0-from lesser error, 1-fixed limit, 2-207Pb/206Pb, 3-206Pb/238U)
        #best age fixed limit cutoff
        g_matrix_inputs = [
            [4, 10, 0, 5, 0, 0, 1000],   #test 0
            [20, 10, 0, 5, 0, 0, 1000],  #test 1
            [1, 1, 0, 5, 0, 0, 1000],    #test 2
            [20, 10, 1, 1.2, 0, 0, 1000],  #test 3
            [20, 10, 1, 2.5, 1, 0, 1000],  # test 4
            [5, 5, 1, 2.5, 0, 1, 333],  # test 5
            [5, 5, 1, 2.5, 0, 1, 1500],  # test 6
            [5, 5, 0, 5, 0, 2, 1000],  # test 7
            [5, 5, 1, 2.5, 0, 3, 1000]  # test 8
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
                [122, 725, 0.41, 97, 14388, 2665, 404], #file 0, test 4
                [131, 2004, 0.78, 105, 4608, 3265, 615],  # file 0, test 5
                [165, 879, 0.39, 115, 22343, 3265, 404],  # file 0, test 6
                [166, 1969, 0.77, 98, 4178, 3265, 305],  # file 0, test 7
                [165, 688, 0.41, 80, 9578, 3533, 404]  # file 0, test 8
            ],
            [#file 1
                [104, 267, 0.4, 31, 1604, 1081, 133], #file 1, test 0
                [146, 258, 0.32, 22, 1289, 2695, 133],#file 1, test 1
                [44, 311, 0.73, 63, 1895, 1081, 177], #file 1, test 2
                [13, 236, 0.85, 69, 1406, 2695, 160], #file 1, test 3
                [15, 260, 0.87, 62, 1113, 767, 160],  #file 1, test 4
                [49, 216, 0.41, 15, 341, 331, 133],  # file 1, test 5
                [115, 267, 0.37, 28, 1480, 1081, 133],  # file 1, test 6
                [115, 623, 4.95, 64, 43, 1105, 55],  # file 1, test 7
                [115, 267, 0.37, 28, 1480, 1081, 133]  # file 1, test 8
            ]
        ]

    def setUp(self):
        self.set_vars()

    def test_matrices_file(self):
        file_name = 'logs/'+str(datetime.datetime.now())+'_gui'
        file_name = file_name.replace(':', '_')
        file_name = file_name + '.log'
        for file_num in range(len(g_matrix_filenames)):
            with open(file_name, 'a') as the_file:
                the_file.write(str(g_matrix_filenames[file_num])+': \n')
            g_operwindow.open_and_load_file(g_matrix_filenames[file_num], False)
            self.set_vars()

            for test_num in range(len(g_matrix_inputs)):
                with open(file_name, 'a') as the_file:
                    the_file.write('Test #'+str(test_num)+' input: ')
                gui_support.onChange(6, g_matrix_inputs[test_num][0], g_pars_onChange)
                gui_support.onChange(7, g_matrix_inputs[test_num][1], g_pars_onChange)
                gui_support.onChange(5, g_matrix_inputs[test_num][2], g_pars_onChange, g_matrix_inputs[test_num][3],
                                     g_operwindow.entErrFilter, g_operwindow.chbInclude207235Err)
                gui_support.onChange(20,  g_matrix_inputs[test_num][3], g_pars_onChange)
                gui_support.onChange(22, g_matrix_inputs[test_num][4], g_pars_onChange)

                gui_support.onChange(3, g_matrix_inputs[test_num][5], g_pars_onChange, g_matrix_inputs[test_num][6],
                                     g_operwindow.entAgeCutoff)
                gui_support.onChange(19, g_matrix_inputs[test_num][6], g_pars_onChange)

                g_operwindow.clear_and_plot()
                self.assertEqual(main_gui.g_number_of_good_grains[0], g_matrix_results[file_num][test_num][0])
                self.assertEqual(int(main_gui.g_number_of_good_grains[1]), g_matrix_results[file_num][test_num][1])
                self.assertEqual(round(main_gui.g_number_of_good_grains[2], 2), g_matrix_results[file_num][test_num][2])
                self.assertEqual(int(main_gui.g_number_of_good_grains[3]), g_matrix_results[file_num][test_num][3])
                self.assertEqual(int(main_gui.g_number_of_good_grains[4]), g_matrix_results[file_num][test_num][4])
                self.assertEqual(main_gui.g_number_of_good_grains[5], g_matrix_results[file_num][test_num][5])
                self.assertEqual(main_gui.g_number_of_good_grains[6], g_matrix_results[file_num][test_num][6])
                with open(file_name, 'a') as the_file:
                    the_file.write(str(g_matrix_inputs[test_num]) + '\n')
                    the_file.write('Output: ' + str(g_matrix_results[file_num][test_num]) + '\n')
        with open(file_name, 'a') as the_file:
                    the_file.write('Success!')

if __name__ == '__main__':
    unittest.main()