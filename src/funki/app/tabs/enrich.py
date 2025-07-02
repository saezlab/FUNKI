import tkinter as tk
from tkinter import ttk
import pandas as pd
import decoupler as dc
import matplotlib.pyplot as plt

from funki.pipelines import enrichment_analysis

from funki.app.utils import LabeledWidget
from funki.app.utils import Figure


class TabEnrich(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=0)
        self.rowconfigure(3, weight=0)
        self.rowconfigure(4, weight=1)

        # Gene set panel
        self.net = pd.DataFrame(columns=['source', 'target', 'weight'])
        ttk.Label(
            self,
            text='Gene set collection:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=0, sticky='NSWE')

        # - Gene set collection selector
        self.gset = tk.StringVar()
        self.gset.set('Hallmark')
        combox_gset = LabeledWidget(
            self,
            ttk.Combobox,
            'Gene set collection:',
            lpos='W',
            wget_kwargs={
                'textvariable': self.gset,
                'values': (
                    ['Hallmark', 'CollecTRI', 'DoRothEA', 'PROGENy']
                    + dc.op.show_resources()['name'].to_list()
                )
            },
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        combox_gset.grid(row=1, column=0, sticky='NSWE')
        combox_gset.wg.bind('<<ComboboxSelected>>', self._get_resource)

        # - Organism selector
        self.org = tk.StringVar()
        self.org.set('Human')
        combox_org = LabeledWidget(
            self,
            ttk.Combobox,
            'Organism:',
            lpos='W',
            wget_kwargs={
                'textvariable': self.org,
                'values': (
                    ['Human']
                    + [i.capitalize() for i in dc.op.show_organisms()]
                )
            },
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        combox_org.grid(row=2, column=0, sticky='NSWE')
        combox_org.wg.bind('<<ComboboxSelected>>', self._get_resource)

        # - View button
        self.button_compute = ttk.Button(
            self,
            text='View',
            command=lambda: self.controller.view_data(dtype='gsc'),
        )
        self.button_compute.grid(row=3, column=0, sticky='W')


        # Enrichment panel
        ttk.Label(
            self,
            text='Enrichment:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=1, sticky='NSWE')

        # - Method selector
        self.method = tk.StringVar()
        self.method.set('ULM (Univariate Linear Model)')
        combox_org = LabeledWidget(
            self,
            ttk.Combobox,
            'Method:',
            lpos='W',
            wget_kwargs={
                'textvariable': self.method,
                'width': 50,
                'values': [
                    'AUCell (Area Under the Curve for set enrichment within '
                        'single cells)',
                    'GSEA (Gene Set Enrichment Analysis)',
                    'GSVA (Gene Set Variation Analysis)',
                    'MDT (Multivariate Decision Trees)',
                    'MLM (Multivariate Linear Model)',
                    'ORA (Over Representation Analysis)',
                    'UDT (Univariate Decision Tree)',
                    'ULM (Univariate Linear Model)',
                    'VIPER (Virtual Inference of Protein-activity by Enriched '
                        'Regulon analysis)',
                    'WAGGR (Weighted Aggregate)',
                    'Z-score',
                ]
            },
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        combox_org.grid(row=1, column=1, sticky='NSWE')

        # - Enrichment variable selector
        self.contrast = tk.StringVar()
        self.combox_obs = LabeledWidget(
            self,
            ttk.Combobox,
            'Choose contrast to enrich from: ',
            lpos='w',
            wget_kwargs={
                'state': 'disabled',
                'textvariable': self.contrast,
            },
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_obs.grid(row=2, column=1, sticky='NSEW')

        # - Compute button
        self.button_compute = ttk.Button(
            self,
            text='Enrich',
            command=self.compute,
            state='disabled',
        )
        self.button_compute.grid(row=3, column=1, sticky='W')

        # Figure
        self.fig, self.ax = plt.subplots()

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=4, columnspan=2, sticky='NSWE')


    def _update(self, *ev):

        # Gene set collection
        if self.net.empty:

            self._get_resource()

        # Enrichment variable
        if self.controller.data:

            if 'diff_exp' in self.controller.data.uns['funki']:

                contrasts = sorted(
                    self.controller.data.uns['funki']['diff_exp'].keys()
                )

                if contrasts:

                    contrast = self.contrast.get() or contrasts[0]

                    self.combox_obs.wg.configure(
                            state='readonly',
                            values=contrasts,
                        )
                    self.contrast.set(contrast)

                    # Activate compute button
                    self.button_compute.configure(state='normal')


    def _get_resource(self, *ev):

        res = self.gset.get()
        org = self.org.get().lower()

        # Decoupler has specific method for resource
        if hasattr(dc.op, res.lower()):

            self.net = getattr(dc.op, res.lower())(organism=org)

        else:

            self.net = dc.op.resource(res, organism=org)


    def compute(self, *ev):

        method = self.method.get().split(' (')[0].replace('-', '').lower()

        self.fig, self.ax = plt.subplots()
        self.figframe.newfig(self.fig)

        enrichment_analysis(
            self.controller.data,
            self.net,
            contrast=self.contrast.get(),
            method=method,
            ax=self.ax,
            gset=self.gset.get(),
            org=self.org.get(),
        )
        self.fig.tight_layout()
        self.figframe._update()
