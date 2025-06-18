from tkinter import ttk

from funki.common import _colors


def load_style(root):

    style = ttk.Style(root)
    # Global
    style.configure(
        '.',
        font=('Arial', 12),
        background=_colors['lightgray'],
    )
    # Label
    style.configure(
        'Title.TLabel',
        font=('Arial', 18, 'bold'),
        foreground=_colors['blue'],
    )
    # Table (Treeview)
    style.configure(
        'Treeview',
        font=('Courier', 11),
    )
    # Tabs
    style.configure(
        'TNotebook',
        tabposition='wn',
    )
    style.configure(
        'TNotebook.Tab',
        padding=(5, 5, 5, 5),
        background=_colors['blue'],
        foreground=_colors['white'],
        width=12,
        anchor='EW',
    )
    style.map(
        'TNotebook.Tab',
        background=[
            ('selected', _colors['teal']),
            ('active', _colors['aqua']),
            #('disabled', _colors['teal'])
        ],
        foreground=[
            ('selected', _colors['yellow']),
            ('active', _colors['white']),
            #('disabled', _colors['yellow'])
        ],
    )
    # Combobox
    style.map(
        'TCombobox',
        fieldbackground=[
            ('disabled', _colors['lightgray']),
            ('readonly', _colors['white']),
        ],
        selectbackground=[
            ('readonly', _colors['aqua'])
        ],
    )

    return style
