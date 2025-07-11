import re
import traceback
from threading import Thread

import tkinter as tk
from tkinter import ttk
import pandas as pd
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

from funki.common import _colors


PATH_LOGO = 'src/funki/app/assets/funki_logo.png'
PATH_MSG = 'src/funki/app/assets/msg_home.txt'
PATH_ICON_SWP = 'src/funki/app/assets/icon_swap.png'


def read_text(path):

    with open(path, 'r') as f:

        txt = ''.join(f.readlines())

    return txt


def check_num(n):

    return re.match('^[0-9.]*$', n) is not None


class ProgressBar(ttk.Progressbar):

    def grid(self, **kwargs):

        self.grid_kwargs = kwargs


    def __enter__(self):

        super().grid(**self.grid_kwargs)
        self['value'] = 0
        self.start()
        self.running = True
        self.thread = Thread(target=self._update)
        self.thread.start()
        self.update()


    def __exit__(self, type, value, traceback):

        self.stop()
        self.running = False
        self.grid_forget()


    def _update(self):

        while self.running:

            self.step()


class LabeledWidget(ttk.Frame):
    '''
    Meta-widget consisting of a frame with a label in a specified position
    relative to a desired widget.
    '''

    def __init__(
        self,
        parent,
        Widget,
        label_text,
        lpos='N',
        wget_kwargs={},
        label_kwargs={},
        wget_grid_kwargs={},
        label_grid_kwargs={},
        **options
    ):

        super().__init__(parent, **options)

        self.wg = Widget(self, **wget_kwargs)
        self.label = ttk.Label(self, text=label_text, **label_kwargs)

        # Positioning
        lpos = lpos.lower()[0]

        if lpos not in 'nswe':

            raise ValueError(
                'Label position invalid, must be one of N, S, W or E.'
            )

        wwg = wget_grid_kwargs.pop('weight', 1)
        wlabel = label_grid_kwargs.pop('weight', 0)

        if lpos in 'ns':

            self.rowconfigure(0, weight=wwg if lpos == 's' else wlabel)
            self.rowconfigure(1, weight=wwg if lpos == 'n' else wlabel)

            self.wg.grid(column=0, row=int(lpos == 'n'), **wget_grid_kwargs)
            self.label.grid(column=0, row=int(lpos == 's'), **label_grid_kwargs)

        elif lpos in 'we':

            self.columnconfigure(0, weight=wwg if lpos == 'e' else wlabel)
            self.columnconfigure(1, weight=wwg if lpos == 'w' else wlabel)

            self.wg.grid(row=0, column=int(lpos == 'w'), **wget_grid_kwargs)
            self.label.grid(row=0, column=int(lpos == 'e'), **label_grid_kwargs)


class Figure(ttk.Frame):
    '''
    Container for a matplotlib figure with corresponding toolbar.
    '''

    def __init__(self, parent, fig, **options):

        super().__init__(parent, relief='sunken', borderwidth=1, **options)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)

        self.canvas = FigureCanvasTkAgg(fig, master=self)

        toolbar = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        toolbar.update()
        toolbar.grid(column=0, row=1, sticky='NSWE')

        self.canvas.get_tk_widget().grid(column=0, row=0, sticky='NSWE')


    def _update(self):

        self.canvas.draw()


    def newfig(self, fig):

        self.canvas = FigureCanvasTkAgg(fig, master=self)

        toolbar = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        toolbar.update()
        toolbar.grid(column=0, row=1, sticky='NSWE')

        self.canvas.get_tk_widget().grid(column=0, row=0, sticky='NSWE')


class Table(ttk.Frame):
    '''
    Generates a ttk.Frame displaying a table with values extracted from a pandas
    DataFrame
    '''

    maxcols = 100
    maxrows = 100
    decimals = 3
    charlim = 25

    def __init__(self, parent, df=None, **options):

        super().__init__(parent, **options)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)

        # Setting up
        self.tree = ttk.Treeview(self)
        self.tree.grid(column=0, row=0, sticky='NSEW')

        xscrollbar = ttk.Scrollbar(
            self,
            orient='horizontal',
            command=self.tree.xview,
        )
        xscrollbar.grid(column=0, row=1, sticky='EW')
        yscrollbar = ttk.Scrollbar(
            self,
            orient='vertical',
            command=self.tree.yview,
        )
        yscrollbar.grid(column=1, row=0, sticky='NS')

        self.tree.configure(
            xscrollcommand=xscrollbar.set,
            yscrollcommand=yscrollbar.set,
        )

        if df is not None:

            self.populate(df)


    def populate(self, df):
        '''
        Populates the table with the provided data.
        '''

        if not isinstance(df, pd.DataFrame):

            raise TypeError('Data provided is not a `pandas.DataFrame`')

        cropdf = self._crop_df(df)

        # Setting up contents
        columns = cropdf.columns.astype(str).to_list()
        self.tree['columns'] = columns
        self.tree.column('#0', width=150, anchor='w', stretch=False)

        for c in columns:

            self.tree.heading(c, text=c)
            self.tree.column(c, width=50, anchor='e', stretch=False)

        for i, (n, r) in enumerate(cropdf.iterrows()):

            self.tree.insert('', 'end', text=n, tag=i % 2, values=[
                self._fmt_val(v) for v in r
            ])

        self.tree.tag_configure(0, background=_colors['white'])
        self.tree.tag_configure(1, background=_colors['lightgray'])


    def _crop_df(self, df):
        '''
        Crops the data table based on maximum number of columns/rows for preview
        purposes.
        '''

        # Cropping rows
        if df.shape[0] > self.maxrows:

            top = bttm = self.maxrows // 2
            top += self.maxrows % 2

            x = pd.DataFrame(
                [['...'] * df.shape[1]],
                columns=df.columns,
                index=['...'],
            )

            df = pd.concat([df.iloc[:top, :], x, df.iloc[-bttm:, :]], axis=0)

        # Cropping columns
        if df.shape[1] > self.maxcols:

            left = right = self.maxcols // 2
            left += self.maxcols % 2

            x = pd.DataFrame(
                ['...'] * df.shape[0],
                columns=['...'],
                index=df.index,
            )

            df = pd.concat([df.iloc[:, :left], x, df.iloc[:, -right:]], axis=1)

        return df


    def _fmt_val(self, val):
        '''
        Formats a value in the table based on character/decimal limits. If the
        string is over the limit, cuts the characters to the limit and adds
        ellipsis. If the number has more decimals than the limit, it's rounded
        and cut up to those decimal positions.
        '''

        if isinstance(val, float):

            val = np.round(val, decimals=self.decimals)

        val = str(val)

        if len(val) > self.charlim:

            val = val[:self.charlim - 3] + '...'

        return val


class PopUpTable(tk.Toplevel):
    '''
    Pops up a new window displaying a given table.
    '''

    def __init__(self, df, title=None, save_command=None):

        super().__init__()

        self.title(title if title else '')
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        frame = ttk.Frame(self, padding=(5, 5, 5, 5))
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(0, weight=1)
        frame.rowconfigure(1, weight=0)
        frame.grid(sticky='NSWE')

        if save_command:

            button_save = ttk.Button(
                frame,
                text='Save',
                command=save_command,
                padding=(5, 5, 5, 5),
                state='disabled',
            )
            button_save.grid(row=1, column=0, sticky='W')

        button_close = ttk.Button(
            frame,
            text='Close',
            command=self.destroy,
            padding=(5, 5, 5, 5),
        )
        button_close.grid(row=1, column=1, sticky='E')

        if not df.empty:

            table = Table(frame, df, padding=(5, 5, 5, 5))
            table.grid(row=0, columnspan=2, sticky='NSWE')

            if save_command:

                button_save.configure(state='normal')


class PopUpError(tk.Toplevel):

    def __init__(self, exc, val, tb):

        super().__init__()

        name = exc.__name__
        self.title(f'FUNKI - {name}!')

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        frame = ttk.Frame(self, padding=(5, 5, 5, 5))

        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        frame.rowconfigure(1, weight=0)

        frame.grid(row=0, column=0, sticky='NSEW')

        text = ScrollText(
            frame,
            text='\n'.join(traceback.format_tb(tb)) + '\n' + f'{name}: {val}'
        )
        text.grid(column=0, row=0, sticky='NSEW')

        button_close = ttk.Button(
            frame,
            text='Close',
            command=self.destroy,
            padding=(5, 5, 5, 5),
        )
        button_close.grid(row=1, column=0)


class ScrollText(ttk.Frame):

    def __init__(self, parent, text=None, **options):

        super().__init__(parent)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        #self.rowconfigure(1, weight=0)

        self.text = tk.Text(self, **options)
        self._update(text)
        self.text.grid(row=0, column=0, sticky='NSWE')

        #xscrollbar = ttk.Scrollbar(
        #    self,
        #    orient='horizontal',
        #    command=self.text.xview,
        #)
        #xscrollbar.grid(column=0, row=1, sticky='EW')
        yscrollbar = ttk.Scrollbar(
            self,
            orient='vertical',
            command=self.text.yview,
        )
        yscrollbar.grid(column=1, row=0, sticky='NS')

        self.text.configure(
            #xscrollcommand=xscrollbar.set,
            yscrollcommand=yscrollbar.set,
        )


    def _update(self, text):

        if text:

            self.text.configure(state='normal')

            self.text.delete('0.0', 'end')
            self.text.insert('end', text)
            self.text.configure(state='disabled')
