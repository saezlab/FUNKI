import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt

from funki.pipelines import sc_quality_control

from utils import Figure

class TabNorm(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=0)
        self.rowconfigure(3, weight=0)
        self.rowconfigure(4, weight=0)
        self.rowconfigure(5, weight=1)

        self.fig_raw, self.ax_raw = plt.subplots(nrows=2, ncols=3)

        self.figframe_raw = Figure(self, self.fig_raw)
        self.figframe_raw.grid(row=5, columnspan=2, sticky='NSWE')

    
    def update(self, *ev):

        if self.controller.data:

            for ax in self.ax_raw.flat:

                ax.clear()

            sc_quality_control(self.controller.data, ax=self.ax_raw)
            self.figframe_raw.update()