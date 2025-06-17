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

        # Raw data panel
        title_raw = ttk.Label(self, text='Raw data:', style='Title.TLabel')
        title_raw.grid(row=0, column=0, sticky='NSWE')

        # Obs data panel
        title_obs = ttk.Label(self, text='Metadata:', style='Title.TLabel')
        title_obs.grid(row=0, column=1, sticky='NSWE')
        self.combox = LabeledWidget(
            self,
            ttk.Combobox,
            'Select variable to visualize: ',
            lpos='w',
            wget_kwargs={'state': 'disabled'}
        )
        self.combox.grid(row=1, column=1)
        self.combox.wg.bind('<<ComboboxSelected>>', self._update)

        self.fig, self.ax = plt.subplots()

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=2, column=1, sticky='NSWE')


    def _update(self, *ev):

        if self.controller.data:

            if not self.controller.data.obs.empty:

                # Set combox
                if not self.combox.wg.get():

                    obs_keys = list(self.controller.data.obs_keys())

                    self.combox.wg.configure(
                        state='readonly',
                        values=obs_keys,
                    )
                    self.combox.wg.set(obs_keys[0])

                # Plot
                self.ax.clear()
                plot_obs(
                    self.controller.data,
                    obs_var=self.combox.wg.get(),
                    ax=self.ax
                )
                self.figframe._update()
