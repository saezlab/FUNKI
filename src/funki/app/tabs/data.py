import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from funki.input import DataSet
from funki.plots import plot_obs
from funki.preprocessing import sc_pseudobulk

from utils import Figure
from utils import LabeledWidget


class TabData(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=0)
        self.rowconfigure(3, weight=0)
        self.rowconfigure(4, weight=0)
        self.rowconfigure(5, weight=0)
        self.rowconfigure(6, weight=1)

        # Raw data panel
        ttk.Label(
            self,
            text='Raw data:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=0, sticky='NW')

        self.nobs = tk.StringVar()
        self.nvar = tk.StringVar()

        LabeledWidget(
            self,
            ttk.Label,
            'No. of observations: ',
            lpos='w',
            wget_kwargs={'textvariable': self.nobs},
            wget_grid_kwargs={'sticky': 'w', 'weight': 0},
            label_grid_kwargs={'sticky': 'w', 'weight': 0},
        ).grid(row=1, column=0, sticky='NW')

        LabeledWidget(
            self,
            ttk.Label,
            'No. of variables: ',
            lpos='w',
            wget_kwargs={'textvariable': self.nvar},
            wget_grid_kwargs={'sticky': 'w', 'weight': 0},
            label_grid_kwargs={'sticky': 'w', 'weight': 0},
        ).grid(row=2, column=0, sticky='NW')

        # Pseudobulking panel
        ttk.Label(
            self,
            text='Pseudobulking:',
            style='Title.TLabel',
            width=20
        ).grid(row=3, column=0, sticky='NW')

        self.sample = tk.StringVar()
        self.combox_pb_sample = LabeledWidget(
            self,
            ttk.Combobox,
            'Sample variable to pseudobulk:',
            lpos='n',
            wget_kwargs={'state': 'disabled', 'textvariable': self.sample},
            wget_grid_kwargs={'sticky': 'WE', 'weight': 1},
            label_grid_kwargs={'sticky': 'W', 'weight': 1},
        )
        self.combox_pb_sample.grid(row=4, column=0, sticky='NW')
        self.combox_pb_sample.wg.bind('<<ComboboxSelected>>', self._update)

        self.group = tk.StringVar()
        self.combox_pb_group = LabeledWidget(
            self,
            ttk.Combobox,
            'Grouping variable for pseudobulk:',
            lpos='n',
            wget_kwargs={'state': 'disabled', 'textvariable': self.group},
            wget_grid_kwargs={'sticky': 'WE', 'weight': 1},
            label_grid_kwargs={'sticky': 'W', 'weight': 1},
        )
        self.combox_pb_group.grid(row=5, column=0, sticky='NW')
        self.combox_pb_group.wg.bind('<<ComboboxSelected>>', self._update)

        # Compute button
        self.button_compute = ttk.Button(
            self,
            text='Apply',
            command=self.compute,
            state='disabled',
        )
        self.button_compute.grid(row=6, column=0, sticky='NW')

        # Obs data panel
        ttk.Label(
            self,
            text='Metadata:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=1, sticky='NSWE')

        self.obs_key = tk.StringVar()
        self.combox_obs = LabeledWidget(
            self,
            ttk.Combobox,
            'Select variable to visualize: ',
            lpos='w',
            wget_kwargs={'state': 'disabled', 'textvariable': self.obs_key},
            wget_grid_kwargs={'sticky': 'WE', 'weight': 0},
            label_grid_kwargs={'sticky': 'W', 'weight': 0},
        )
        self.combox_obs.grid(row=1, column=1, sticky='NSWE')
        self.combox_obs.wg.bind('<<ComboboxSelected>>', self._update)

        self.fig, self.ax = plt.subplots()

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=2, column=1, rowspan=5, sticky='NSWE')


    def _update(self, *ev):

        if self.controller.data:

            no, nv = self.controller.data.shape
            self.nobs.set(str(no))
            self.nvar.set(str(nv))

            if not self.controller.data.obs.empty:

                # Comboboxes for pseudobulk
                obs_keys_str = sorted([
                    c for c in self.controller.data.obs_keys()
                    if all([
                        isinstance(i, str)
                        for i in self.controller.data.obs[c]
                    ])
                ])

                if obs_keys_str:

                    self.combox_pb_sample.wg.configure(
                            state='readonly',
                            values=obs_keys_str,
                    )
                    self.sample.set(self.sample.get() or obs_keys_str[0])
                    self.combox_pb_group.wg.configure(
                            state='readonly',
                            values=obs_keys_str,
                    )
                    self.group.set(self.group.get() or obs_keys_str[0])

                    self.button_compute.configure(state='normal')

                # Set combobox for viz
                obs_keys = sorted(self.controller.data.obs_keys())

                if obs_keys:

                    self.combox_obs.wg.configure(
                        state='readonly',
                        values=obs_keys,
                    )
                    self.obs_key.set(self.obs_key.get() or obs_keys[0])

                    # Plot
                    self.ax.clear()
                    plot_obs(
                        self.controller.data,
                        obs_var=self.obs_key.get(),
                        ax=self.ax
                    )
                    self.figframe._update()


    def compute(self, *ev):

        aux = self.controller.data.copy()
        sc_pseudobulk(
            self.controller.data,
            self.sample.get(),
            groups_col=self.group.get()
        )
        # TODO: Should make this a standard function
        self.controller.data = DataSet(
            self.controller.data.uns['pseudobulk'].copy()
        )
        self.controller.data.uns['nonpseudo_original'] = aux.copy()

        self.controller._update()
