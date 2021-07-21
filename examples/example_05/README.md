# Running `calphy` from jupyter notebooks

In this example, `calphy` will be used as a library to run Example 01 directly from a jupyter notebook. Please check example 01 before completing this example. We start by import a function to read the input file.


```python
from calphy.input import read_inputfile
```


```python
options = read_inputfile("input.yaml")
```

We can check the keys of the options dictionary.


```python
options.keys()
```




    dict_keys(['calculations', 'md', 'queue', 'conv', 'element', 'mass', 'nelements'])



The individual methods that are required to run the calculation can be imported from the `queuekernel` module.


```python
import calphy.queuekernel as cq 
```

First, we set up a class which prepares everything for the calculation. It takes `options` as the first argument, followed by an integer argument. The second argument indicates which calculation from the calculation block is run.


```python
job = cq.setup_calculation(options, 0)
```


```python
job
```




    solid system with T=100.000000, P=0.000000 in bcc lattice for mode fe



The specifics of the calculation can be obtained by,


```python
job.calc
```




    {'mode': 'fe',
     'temperature': 100,
     'pressure': 0,
     'lattice': 'BCC',
     'state': 'solid',
     'temperature_stop': 100,
     'nelements': 1,
     'element': ['Fe'],
     'lattice_constant': 2.8842,
     'iso': True,
     'fix_lattice': False,
     'repeat': [5, 5, 5],
     'nsims': 1,
     'thigh': 200.0,
     'directory': 'fe-BCC-100-0'}



These properties can also be accessed individually.


```python
job.t, job.p, job.l
```




    (100, 0, 'bcc')



Now finally the calculation can be run


```python
job = cq.run_calculation(job)
```

The results can be obtained through the `report` variable


```python
job.report
```




    {'input': {'temperature': 100,
      'pressure': 0.0,
      'lattice': 'bcc',
      'element': 'Fe',
      'concentration': '1'},
     'average': {'vol/atom': 11.994752673986264, 'spring_constant': '3.35'},
     'results': {'free_energy': -4.263447380763835,
      'error': 0.0,
      'reference_system': 0.01514568505240216,
      'work': -4.278593065816237,
      'pv': 0.0}}



or individually as class attributes


```python
job.fe, job.w, job.pv
```




    (-4.263447380763835, -4.278593065816237, 0)



If more than one calculation is present, they should be run individually. For example, we use the `input2.yaml` file.


```python
options = read_inputfile("input2.yaml")
```

This input file has two structures: BCC and liquid. We can check the list of calculations


```python
len(options["calculations"])
```




    2



If you want to run the second calculation in the list, we have to set up the job as follows.


```python
job = cq.setup_calculation(options, 1)
```
