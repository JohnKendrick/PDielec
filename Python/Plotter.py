#!/usr/bin/python
"""Do some plotting for the dielectric functions"""
from __future__ import print_function, division
import sys
import numpy as np


class Plotter:
    """Do some plotting for the dielectric functions"""
    def __init__(self):
        self.methods = []
        self.volume_fractions = []
        self.sizes = []
        self.shapes = []
        self.shape_data = []
        self.plot_numbers = []
        self.traces = []
        self.absorption_coefficients = []
        self.molar_absorption_coefficients = []
        self.frequencies = []
        return

    def add_dielectric(self, nplot, method, vf_type, size, shape, data, frequency, trace, absorption, molar_absorption):
        """Add the trace of a dielectric to be plotted, add the absorption and the molar absorption"""
        if nplot not in self.plot_numbers:
            self.plot_numbers.append(nplot)
            self.methods.append(method)
            self.volume_fractions.append(vf_type)
            self.sizes.append(size)
            self.shapes.append(shape)
            self.shape_data.append(data)
            self.traces.append([])
            self.absorption_coefficients.append([])
            self.molar_absorption_coefficients.append([])
            self.frequencies.append([])
        self.frequencies[nplot].append(frequency)
        self.traces[nplot].append(trace)
        self.absorption_coefficients[nplot].append(absorption)
        self.molar_absorption_coefficients[nplot].append(molar_absorption)
        return

    def printout(self, fd_csvfile):
        """Print data to a csv file"""
        if fd_csvfile == 0:
            return
        # Print the header
        output = ",,,,,cm-1,thz"
        for method, shape, data, vf_type, size in zip(self.methods, self.shapes, self.shape_data, self.volume_fractions, self.sizes):
            strdata = str(data)
            strdata = strdata.replace(" ", "")
            strdata = strdata.replace(",", " ")
            strdata = strdata.replace("'", "")
            strdata = strdata.replace('"', "")
            ssize=""
            if size > 1.0e-6:
                ssize = "/{:f}".format(size).rstrip("0").rstrip(".")+"mu"
            output += ","+method+"_"+shape+strdata+"_real_dielectric"+"("+vf_type+ssize+")"
            output += ","+method+"_"+shape+strdata+"_imag_dielectric"+"("+vf_type+ssize+")"
            output += ","+method+"_"+shape+strdata+"_absorption_coeficient_cm-1"+"("+vf_type+ssize+")"
            output += ","+method+"_"+shape+strdata+"_molar_absorption_coeficient_cm-1.L.mole-1"+"("+vf_type+ssize+")"
            # end loop over method.....
        print(output, file=fd_csvfile)
        # Use the first frequency list
        for iv, v in enumerate(self.frequencies[0]):
            thz = v * 0.0299792458
            output = ",,,,,{:12f},{:12f}".format(v, thz)
            for dri, absorption, molar_absorption in zip(self.traces, self.absorption_coefficients, self.molar_absorption_coefficients):
                trace_real = float(np.real(dri[iv]))
                trace_imag = float(np.imag(dri[iv]))
                absrp = absorption[iv]
                mabsrp = molar_absorption[iv]
                output += "," + ",".join("{:20.8f}".format(p) for p in [trace_real, trace_imag, absrp, mabsrp])
                # end loop dri, absorption
            print(output, file=fd_csvfile)
        # end loop over frequency
        return

    def plot(self, types):
        """Plot driver routine"""
        for plot in types:
            if plot == 'imaginary':
                self.plot_imaginary()
            elif plot == 'real':
                self.plot_real()
            elif plot == 'absorption':
                self.plot_absorption()
            elif plot == 'molar_absorption':
                self.plot_molar_absorption()
            elif plot == 'molarAbsorption':
                self.plot_molar_absorption()
            elif plot == 'molarabsorption':
                self.plot_molar_absorption()
            else:
                print('Unkown plot type {}'.format(plot))
                exit(1)
            # endif
        # endfor
        return

    def plot_molar_absorption(self):
        """Plot loss data """
        import pylab as pl
        # Choose a suitable scale TO has to be treated differently, TO was added last
        maxd = 0.0
        max_to = 1.2*np.amax(self.molar_absorption_coefficients[-1])
        # loop over all methods except the last, TO
        for molar_absorption in self.molar_absorption_coefficients[:-1]:
            maxd = max(maxd, np.amax(molar_absorption))
        # end loop over entries
        to_scale = maxd/max_to
        to_scale = 1.0
        molar_absorption = np.array(self.molar_absorption_coefficients[-1])
        self.molar_absorption_coefficients[-1] = to_scale * molar_absorption
        plot_list = []
        plot_labels = []
        for method, shape, data, frequencies, molar_absorption, vf_type, size in zip(self.methods, self.shapes, self.shape_data, self.frequencies, self.molar_absorption_coefficients, self.volume_fractions, self.sizes):
            ssize=""
            if size > 1.0e-6:
                ssize = "/{:f}".format(size).rstrip("0").rstrip(".")+"mu"
            label = method+"_"+shape+str(data)+"("+vf_type+ssize+")"
            label = label.replace(" ", "")
            label = label.replace(",", " ")
            label = label.replace("'", "")
            label = label.replace('"', '')
            plot, = pl.plot(frequencies, molar_absorption, label=label)
            plot_list.append(plot)
            plot_labels.append(label)
        # end loop over method
        pl.legend(plot_list,  plot_labels, loc='best', numpoints=1)
        pl.xlabel('Frequency (cm-1)')
        pl.ylabel('Molar Absorption (cm-1.L.mole-1)')
        pl.title('Molar Absorption')
        pl.show()
        return

    def plot_absorption(self):
        """Plot absorption data """
        import pylab as pl
        # Choose a suitable scale TO has to be treated differently, TO was added last
        maxd = 0.0
        max_to = 1.2*np.amax(self.absorption_coefficients[-1])
        # loop over all methods except the last, TO
        for absorptions in self.absorption_coefficients[:-1]:
            maxd = max(maxd, np.amax(absorptions))
        # end loop over entries
        to_scale = maxd/max_to
        to_scale = 1.0
        absorptions = np.array(self.absorption_coefficients[-1])
        self.absorption_coefficients[-1] = to_scale * absorptions
        plot_list = []
        plot_labels = []
        for method, shape, data, frequencies, absorptions, vf_type, size in zip(self.methods, self.shapes, self.shape_data, self.frequencies, self.absorption_coefficients, self.volume_fractions, self.sizes):
            ssize=""
            if size > 1.0e-6:
                ssize = "/{:f}".format(size).rstrip("0").rstrip(".")+"mu"
            label = method+"_"+shape+str(data)+"("+vf_type+ssize+")"
            label = label.replace(" ", "")
            label = label.replace(",", " ")
            label = label.replace("'", "")
            label = label.replace('"', '')
            plot, = pl.plot(frequencies, absorptions, label=label)
            plot_list.append(plot)
            plot_labels.append(label)
        # end loop over method
        pl.legend(plot_list,  plot_labels, loc='best', numpoints=1)
        pl.xlabel('Frequency (cm-1)')
        pl.ylabel('Absorption Coefficient (cm-1)')
        pl.title('Absorption Coefficient')
        pl.show()
        return

    def plot_imaginary(self):
        """Plot loss data """
        import pylab as pl
        # Choose a suitable scale TO has to be treated differently, TO was added last
        maxd = 0.0
        max_to = 1.2*np.amax(np.imag(self.traces[-1]))
        # loop over all methods except the last, TO
        for traces in self.traces[:-1]:
            real_trace = np.imag(traces)
            maxd = max(maxd, np.amax(real_trace))
        # end loop over entries
        to_scale = maxd/max_to
        to_scale = 1.0
        traces = np.array(self.traces[-1])
        self.traces[-1] = to_scale * traces
        plot_list = []
        plot_labels = []
        for method, shape, data, frequencies, traces, vf_type, size in zip(self.methods, self.shapes, self.shape_data, self.frequencies, self.traces, self.volume_fractions, self.sizes):
            ssize=""
            if size > 1.0e-6:
                ssize = "/{:f}".format(size).rstrip("0").rstrip(".")+"mu"
            label = method+"_"+shape+str(data)+"("+vf_type+ssize+")"
            label = label.replace(" ", "")
            label = label.replace(",", " ")
            label = label.replace("'", "")
            label = label.replace('"', '')
            plot, = pl.plot(frequencies, np.imag(traces), label=label)
            plot_list.append(plot)
            plot_labels.append(label)
        # end loop over method
        pl.legend(plot_list,  plot_labels, loc='best', numpoints=1)
        pl.xlabel('Frequency (cm-1)')
        pl.ylabel('Imaginary Permittivity')
        pl.title('Imaginary Permittivity')
        pl.show()
        return

    def plot_real(self):
        """Plot real data """
        import pylab as pl
        # Choose a suitable scale TO has to be treated differently, TO was added last
        maxd = 0.0
        max_to = 1.2*np.amax(np.real(self.traces[-1]))
        # loop over all methods except the last, TO
        for traces in self.traces[:-1]:
            real_trace = np.real(traces)
            maxd = max(maxd, np.amax(real_trace))
        # end loop over entries
        to_scale = maxd/max_to
        to_scale = 1.0
        traces = np.array(self.traces[-1])
        self.traces[-1] = to_scale * traces
        plot_list = []
        plot_labels = []
        for method, shape, data, frequencies, traces, vf_type, size in zip(self.methods, self.shapes, self.shape_data, self.frequencies, self.traces, self.volume_fractions, self.sizes):
            ssize=""
            if size > 1.0e-6:
                ssize = "/{:f}".format(size).rstrip("0").rstrip(".")+"mu"
            label = method+"_"+shape+str(data)+"("+vf_type+ssize+")"
            label = label.replace(" ", "")
            label = label.replace(",", " ")
            label = label.replace("'", "")
            label = label.replace('"', '')
            plot, = pl.plot(frequencies, np.real(traces), label=label)
            plot_list.append(plot)
            plot_labels.append(label)
        # end loop over method
        pl.legend(plot_list,  plot_labels, loc='best', numpoints=1)
        pl.xlabel('Frequency (cm-1)')
        pl.ylabel('Real Permittivity')
        pl.title('Real Permittivity')
        pl.show()
        return


def print_ints(title, ints, no_per_line=8, format="{:9d}", file=sys.stdout, separator=" "):
    """Print ints data """
    #
    # Print out a list of ints prettily
    #
    len_ints = len(ints)
    print(" ", file=file)
    if title != "":
        print(title, file=file)
    nlines = int((len_ints - 1) / no_per_line) + 1
    start = 0
    for i in range(nlines):
        end = start + no_per_line
        if end > len_ints:
            end = len_ints
        print(" " + separator.join(format.format(r) for r in ints[start:end]), file=file)
        start = start + no_per_line
    # end for i
    return

def print_strings(title, strings, no_per_line=8, format="{:9s}", file=sys.stdout, separator=" "):
    """Print strings data """
    #
    # Print out a list of strings prettily
    #
    len_strings = len(strings)
    print(" ", file=file)
    if title != "":
        print(title, file=file)
    nlines = int((len_strings - 1) / no_per_line) + 1
    start = 0
    for i in range(nlines):
        end = start + no_per_line
        if end > len_strings:
            end = len_strings
        print(" " + separator.join(format.format(r) for r in strings[start:end]), file=file)
        start = start + no_per_line
    # end for i
    return


def print_reals(title, reals, no_per_line=8, format="{:9.2f}", file=sys.stdout, separator=" "):
    """Print reals data """
    #
    # Print out a list of reals prettily
    #
    len_reals = len(reals)
    print(" ", file=file)
    if title != "":
        print(title, file=file)
    nlines = int((len_reals - 1) / no_per_line) + 1
    start = 0
    for i in range(nlines):
        end = start + no_per_line
        if end > len_reals:
            end = len_reals
        print(" " + separator.join(format.format(r) for r in reals[start:end]), file=file)
        start = start + no_per_line
    # end for i
    return


def print3x3(title, array, format="{:14.6f}", file=sys.stdout, separator=" "):
    """Print 3x3 matrix"""
    #
    # Print out a 3x3 tensor matrix
    #
    print(" ", file=file)
    if title != "":
        print(title, file=file)
    for i in range(3):
        print("      "+" ".join(format.format(p) for p in array[i]), file=file)
    # end for i
    return
