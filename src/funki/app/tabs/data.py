import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt

from funki.pipelines import sc_quality_control

from utils import Figure


class TabData(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.fig, self.ax = plt.subplots(nrows=2, ncols=3)

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=0, column=0, sticky='NSWE')


    def update(self):

        if self.controller.data:

            sc_quality_control(self.controller.data, ax=self.ax)
            self.figframe.update()
