import param
import pandas as pd
import datetime as dt
from standard_workflows import analysis_params

class BaseClass(param.Parameterized):
    x                       = param.Parameter(default=3.14, doc="X position")
    y                       = param.Parameter(default="Not editable", constant=True)
    string_value            = param.String(default="str", doc="A string")
    num_int                 = param.Integer(50000, bounds=(-200, 100000))
    unbounded_int           = param.Integer(23)
    float_with_hard_bounds  = param.Number(8.2, bounds=(7.5, 10))
    float_with_soft_bounds  = param.Number(0.5, bounds=(0, None), softbounds=(0,2))
    unbounded_float         = param.Number(30.01, precedence=0)
    hidden_parameter        = param.Number(2.718, precedence=-1)
    integer_range           = param.Range(default=(3, 7), bounds=(0, 10))
    float_range             = param.Range(default=(0, 1.57), bounds=(0, 3.145))
    dictionary              = param.Dict(default=analysis_params.analysis_params['default']['basepaths'])
    


class Example(BaseClass):
    """An example Parameterized class"""
    timestamps = []

    boolean                 = param.Boolean(True, doc="A sample Boolean parameter")
    color                   = param.Color(default='#FFFFFF')
    date                    = param.Date(dt.datetime(2017, 1, 1),
                                         bounds=(dt.datetime(2017, 1, 1), dt.datetime(2017, 2, 1)))
    dataframe               = param.DataFrame(pd._testing.makeDataFrame().iloc[:3])
    select_string           = param.ObjectSelector(default="yellow", objects=["red", "yellow", "green"])
    select_fn               = param.ObjectSelector(default=list,objects=[list, set, dict])
    int_list                = param.ListSelector(default=[3, 5], objects=[1, 3, 5, 7, 9], precedence=0.5)
    single_file             = param.FileSelector(path='../../*/*.py*', precedence=0.5)
    multiple_files          = param.MultiFileSelector(path='../../*/*.py?', precedence=0.5)
    
    def update_timestamp(x): 
        x.timestamps.append(dt.datetime.utcnow())

    def update_string_value(self, event): 
        self.string_value = str(len(self.timestamps))

    record_timestamp        = param.Action(update_timestamp, doc="""Record timestamp.""", precedence=0.7)




import panel as pn
pn.extension()

base = BaseClass()
example=Example()

pn.Row(example.param, base.param).servable()


import pandas as pd

df = pd.DataFrame({
  'A': [1, 2, 3, 4],
  'B': [10, 20, 30, 40]
})

df_pane = pn.panel(df)

pn.Row(df_pane).servable()


x = pn.widgets.IntSlider(name='x', start=0, end=100)

def square(x):
    example.update_string_value(str(x))
    return f'{x} squared is {x**2} and {example.string_value}'

pn.Row(pn.bind(square, x)).servable()
pn.panel(x).servable()

#pn.Row(pn.bind(example.update_string_value, x)).servable()

#pn.Row(pn.bind(example.string_value, x)).servable()