import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from funki.pipelines import differential_expression

from funki.app.utils import Figure
from funki.app.utils import LabeledWidget
from funki.app.utils import PATH_ICON_SWP


class TabDex(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=0)
        self.rowconfigure(3, weight=0)
        self.rowconfigure(4, weight=0)
        self.rowconfigure(5, weight=1)

        ttk.Label(
            self,
            text='Differential expression:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, columnspan=3, sticky='NSWE')

        # Contrast variable selector
        self.contrast_var = tk.StringVar()
        self.combox_contrast_var = LabeledWidget(
            self,
            ttk.Combobox,
            'Select group variable for the contrast: ',
            lpos='n',
            wget_kwargs={
                'state': 'disabled',
                'textvariable': self.contrast_var,
            },
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_contrast_var.grid(row=1, columnspan=2, sticky='NSEW')
        self.combox_contrast_var.wg.bind('<<ComboboxSelected>>', self._update)

        # TODO: Handle multiple groups?

        # Contrast group A
        self.groupA = tk.StringVar()
        self.combox_A = LabeledWidget(
            self,
            ttk.Combobox,
            'Group of samples to contrast from: ',
            lpos='n',
            wget_kwargs={'state': 'disabled', 'textvariable': self.groupA},
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_A.grid(row=2, column=0, sticky='NSWE')
        self.combox_A.wg.bind('<<ComboboxSelected>>', self._update)

        # Contrast group B
        self.groupB = tk.StringVar()
        self.combox_B = LabeledWidget(
            self,
            ttk.Combobox,
            'Group of samples to contrast against: ',
            lpos='n',
            wget_kwargs={'state': 'disabled', 'textvariable': self.groupB},
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_B.grid(row=2, column=2, sticky='NSWE')
        self.combox_B.wg.bind('<<ComboboxSelected>>', self._update)

        # Swap button
        icon = tk.PhotoImage(file=PATH_ICON_SWP)
        self.button_swap = ttk.Button(
            self,
            image=icon,
            command=self.swap_ab,
            state='disabled',
        )
        self.button_swap.image = icon
        self.button_swap.grid(row=2, column=1, sticky='SW')

        # Thresholds
        self.thr_logfc = tk.DoubleVar(value=1.0)
        LabeledWidget(
            self,
            ttk.Entry,
            'log(FC) threshold: ',
            lpos='w',
            wget_kwargs={
                'textvariable': self.thr_logfc,
                'validate': 'key',
                'validatecommand': self.controller.check_num,
                'width': 5,
            }
        ).grid(column=0, row=3, sticky='W')

        self.thr_pval = tk.DoubleVar(value=0.05)
        LabeledWidget(
            self,
            ttk.Entry,
            'P-value threshold: ',
            lpos='w',
            wget_kwargs={
                'textvariable': self.thr_pval,
                'validate': 'key',
                'validatecommand': self.controller.check_num,
                'width': 5,
            }
        ).grid(column=2, row=3, sticky='W')

        # DEX methods
        method_frame = ttk.Frame(self)
        method_frame.columnconfigure(0, weight=0)
        method_frame.columnconfigure(1, weight=1)
        method_frame.columnconfigure(2, weight=1)

        ttk.Label(
            method_frame,
            text='Select method:'
        ).grid(row=0, column=0, sticky='W')

        self.method = tk.StringVar(value='limma')
        ttk.Radiobutton(
            method_frame,
            text='limma',
            variable=self.method,
            value='limma'
        ).grid(row=0, column=1, sticky='E')
        ttk.Radiobutton(
            method_frame,
            text='PyDESeq2',
            variable=self.method,
            value='pydeseq2'
        ).grid(row=0, column=2, sticky='E')

        method_frame.grid(row=4, columnspan=2, sticky='NSEW')

        # Compute/view button
        frame_buttons = ttk.Frame(self)
        frame_buttons.columnconfigure(0, weight=1)
        frame_buttons.columnconfigure(1, weight=1)

        self.button_compute = ttk.Button(
            frame_buttons,
            text='Compute',
            command=self.compute,
            state='disabled',
        )
        self.button_compute.grid(row=0, column=0, sticky='E')
        self.button_view = ttk.Button(
            frame_buttons,
            text='View',
            command=lambda: self.controller.view_data(dtype='dex'),
            state='disabled',
        )
        self.button_view.grid(row=0, column=1, sticky='E')

        frame_buttons.grid(row=4, column=2, sticky='NSWE')

        # Figure
        self.fig, self.ax = plt.subplots()

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=5, columnspan=3, sticky='NSWE')


    def _update(self, *ev): # TODO: Could probably split to simplify

        obs_key = self.contrast_var.get()

        if self.controller.data and not self.controller.data.obs.empty:

            # Handling available obs_keys
            obs_keys = sorted([
                c for c in self.controller.data.obs_keys()
                if all([
                    isinstance(i, str)
                    for i in self.controller.data.obs[c]
                ])
            ])

            if obs_keys:

                obs_key = obs_key or obs_keys[0]
                self.combox_contrast_var.wg.configure(
                        state='readonly',
                        values=obs_keys,
                    )
                self.contrast_var.set(obs_key)

                # Handling available options in obs_key
                options = sorted(set(self.controller.data.obs[obs_key]))
                curA = self.groupA.get()
                curB = self.groupB.get()

                # A/B are empty or new obs_key
                if (
                    (not curA or not curB)
                    or (curA not in options or curB not in options)
                ):

                    curA = options[0]
                    curB = options[1]

                self.combox_A.wg.configure(
                    state='readonly',
                    values=[n for n in options if n != curB],
                )
                self.groupA.set(curA)

                self.combox_B.wg.configure(
                    state='readonly',
                    values=[n for n in options if n!= curA],
                )
                self.groupB.set(curB)

                # Reactivating buttons
                self.button_compute.configure(state='normal')
                self.button_swap.configure(state='normal')


    def swap_ab(self):

        curA = self.groupA.get()
        curB = self.groupB.get()

        if curA and curB:

            obs_key = self.contrast_var.get()

            options = sorted(set(self.controller.data.obs[obs_key]))

            self.combox_A.wg.configure(
                values=[i for i in options if i != curA],
            )
            self.groupA.set(curB)
            self.combox_B.wg.configure(
                values=[i for i in options if i != curB],
            )
            self.groupB.set(curA)


    def compute(self):

        name = f'{self.groupA.get()}_vs_{self.groupB.get()}'

        self.ax.clear()
        differential_expression(
            self.controller.data,
            self.contrast_var.get(),
            self.groupA.get(),
            self.groupB.get(),
            logfc_thr=self.thr_logfc.get(),
            fdr_thr=self.thr_pval.get(),
            method=self.method.get(),
            ax=self.ax,
        )
        self.fig.tight_layout()
        self.figframe._update()

        # Activating save/view menus and button
        self.button_view.configure(state='normal')
        self.controller.menu_view.entryconfigure(
            'Differential expression...',
            state='normal'
        )
        self.controller.menu_save.entryconfigure(
            'Differential expression...',
            state='normal'
        )

        # Adding submenu for current contrast
        self.controller.menu_view_dex.add_command(
            label=name,
            command=lambda: self.controller.view_data(dtype='dex', key=name)
        )
        self.controller.menu_save_dex.add_command(
            label=name,
            command=lambda: self.controller.save_file(dtype='dex', key=name)
        )

        self.controller._update()
