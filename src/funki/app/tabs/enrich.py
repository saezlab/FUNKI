import tkinter as tk
from tkinter import ttk
import pandas as pd
import decoupler as dc

from funki.analysis import enrich

from utils import LabeledWidget


class TabEnrich(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=1)

        # Gene set panel
        self.net = pd.DataFrame(columns=['source', 'target', 'weight'])
        ttk.Label(
            self,
            text='Gene Set collection:',
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
            }
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
            }
        )
        combox_org.grid(row=2, column=0, sticky='NSWE')
        combox_org.wg.bind('<<ComboboxSelected>>', self._get_resource)

        # Enrichment panel
        ttk.Label(
            self,
            text='Enrichment:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=1, sticky='NSWE')


    def _update(self, *ev):

        self._get_resource()


    def _get_resource(self, *ev):

        res = self.gset.get()
        org = self.org.get().lower()
        
        # Decoupler has specific method for resource
        if hasattr(dc.op, res.lower()):
            
            self.net = getattr(dc.op, res.lower())(organism=org)

        else:

            self.net = dc.op.resource(res, organism=org)
