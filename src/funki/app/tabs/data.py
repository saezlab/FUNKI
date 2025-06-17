import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt

from funki.pipelines import sc_quality_control
from funki.plots import plot_obs

from utils import Figure


class TabData(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=1)

        # Raw data viz

        self.fig_raw, self.ax_raw = plt.subplots(nrows=2, ncols=3)

        self.figframe_raw = Figure(self, self.fig_raw)
        self.figframe_raw.grid(row=1, column=0, sticky='NSWE')

        # Obs data viz

        self.combox_obs = ttk.Combobox(
            self,
            state='disabled'
        )
        self.combox_obs.grid(row=0, column=1)
        self.combox_obs.bind('<<ComboboxSelected>>', self.update)

        self.fig_obs, self.ax_obs = plt.subplots()

        self.figframe_obs = Figure(self, self.fig_obs)
        self.figframe_obs.grid(row=1, column=1, sticky='NSWE')


    def update(self, *ev):

        if self.controller.data:

            for ax in self.ax_raw.flat:

                ax.clear()

            sc_quality_control(self.controller.data, ax=self.ax_raw)
            self.figframe_raw.update()

        if not self.controller.data.obs.empty and self.combox_obs.get():

            self.ax_obs.clear()
            plot_obs(
                self.controller.data,
                obs_var=self.combox_obs.get(),
                ax=self.ax_obs
            )
            self.figframe_obs.update()
