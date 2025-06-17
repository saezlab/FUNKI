import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt

from funki.pipelines import sc_quality_control
from funki.plots import plot_obs

from utils import Figure, LabeledWidget


class TabData(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=1)

        # Raw data viz
        title_raw = ttk.Label(self, text='Raw data QC:', style='Title.TLabel')
        title_raw.grid(row=0, column=0, sticky='NSWE')
        self.fig_raw, self.ax_raw = plt.subplots(nrows=2, ncols=3)

        self.figframe_raw = Figure(self, self.fig_raw)
        self.figframe_raw.grid(row=2, column=0, sticky='NSWE')

        # Obs data viz
        title_obs = ttk.Label(self, text='Metadata:', style='Title.TLabel')
        title_obs.grid(row=0, column=1, sticky='NSWE')
        self.combox_obs = LabeledWidget(
            self,
            ttk.Combobox,
            'Select variable to visualize: ',
            lpos='w',
            wget_kwargs={'state': 'disabled'}
        )
        self.combox_obs.grid(row=1, column=1)
        self.combox_obs.wg.bind('<<ComboboxSelected>>', self.update)

        self.fig_obs, self.ax_obs = plt.subplots()

        self.figframe_obs = Figure(self, self.fig_obs)
        self.figframe_obs.grid(row=2, column=1, sticky='NSWE')


    def update(self, *ev):

        if self.controller.data:

            for ax in self.ax_raw.flat:

                ax.clear()

            sc_quality_control(self.controller.data, ax=self.ax_raw)
            self.figframe_raw.update()

            if not self.controller.data.obs.empty and self.combox_obs.wg.get():

                self.ax_obs.clear()
                plot_obs(
                    self.controller.data,
                    obs_var=self.combox_obs.wg.get(),
                    ax=self.ax_obs
                )
                self.figframe_obs.update()
