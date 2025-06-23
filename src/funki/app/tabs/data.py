import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from funki.plots import plot_obs

from utils import Figure
from utils import LabeledWidget


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
        ttk.Label(
            self,
            text='Raw data:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=0, sticky='NSWE')

        summary = ttk.Frame(self)
        summary.rowconfigure(0, weight=0)
        summary.rowconfigure(1, weight=0)

        self.nobs = tk.StringVar()
        self.nvar = tk.StringVar()

        LabeledWidget(
            summary,
            ttk.Label,
            'No. of observations: ',
            lpos='n',
            wget_kwargs={'textvariable': self.nobs},
            wget_grid_kwargs={'sticky': 'w', 'weight': 0},
            label_grid_kwargs={'sticky': 'w', 'weight': 0},
        ).grid(row=0, sticky='W')

        LabeledWidget(
            summary,
            ttk.Label,
            'No. of variables: ',
            lpos='n',
            wget_kwargs={'textvariable': self.nvar},
            wget_grid_kwargs={'sticky': 'w', 'weight': 0},
            label_grid_kwargs={'sticky': 'w', 'weight': 0},
        ).grid(row=1, sticky='W')

        summary.grid(column=0, row=2, sticky='NSWE')

        # Obs data panel
        ttk.Label(
            self,
            text='Metadata:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=1, sticky='NSWE')

        self.obs_key = tk.StringVar()
        self.combox = LabeledWidget(
            self,
            ttk.Combobox,
            'Select variable to visualize: ',
            lpos='w',
            wget_kwargs={'state': 'disabled', 'textvariable': self.obs_key},
            wget_grid_kwargs={'sticky': 'w', 'weight': 0},
            label_grid_kwargs={'sticky': 'w', 'weight': 0},
        )
        self.combox.grid(row=1, column=1, sticky='NSWE')
        self.combox.wg.bind('<<ComboboxSelected>>', self._update)

        self.fig, self.ax = plt.subplots()

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=2, column=1, sticky='NSWE')


    def _update(self, *ev):

        if self.controller.data:

            no, nv = self.controller.data.shape
            self.nobs.set(str(no))
            self.nvar.set(str(nv))

            if not self.controller.data.obs.empty:

                obs_keys = sorted(self.controller.data.obs_keys())

                # Set combobox
                if obs_keys:

                    obs_key = self.obs_key.get() or obs_keys[0]

                    self.combox.wg.configure(
                        state='readonly',
                        values=obs_keys,
                    )
                    self.obs_key.set(obs_key)

                # Plot
                self.ax.clear()
                plot_obs(
                    self.controller.data,
                    obs_var=obs_key,
                    ax=self.ax
                )
                self.figframe._update()
