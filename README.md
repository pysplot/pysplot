# pysplot: Spatiotemporal Plotting in Python

The `pysplot` package makes spatial visualizations with heliophysics timeseries data.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation

To install the package:
1. Clone package onto local machine and switch to correct branch
```
git clone https://github.com/pysplot/pysplot.git
git checkout whpi-notebook
```
2. Use `pip` to install package:
```
cd /path/to/pysplot
pip install .
```

## Usage

After installation, you can use the package as follows:

```python
from pysplot.io.data import SpatialData, ScienceData, SpatialTimeData
from pysplot.plots.plottypes import orbit_plot
import matplotlib.pyplot as plt

# Initiate data objects
spatial_data = SpatialData(spatial_data_dictionary)
science_data = ScienceData(science_data_dictionary)
combined_splot_data = SpatialTimeData(spatial_data, science=science_data)

# Make plots
fig, ax = plt.subplots()
orbit_plot(combined_splot_data.data, ax, 'x', 'y')
```

For a full example of how to use `pysplot`, refer to our example notebooks.


## Contributing

We welcome contributions to improve this package! To contribute, please follow these steps:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Make your changes.
4. Commit your changes (`git commit -am 'Add new feature'`).
5. Push to the branch (`git push origin feature-branch`).
6. Create a pull request.

Please make sure your code adheres to the existing style guidelines and includes tests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
