"""Different types of spatial plots."""

def orbit_plot(data, ax, x_axis_column, y_axis_column, plotting_kwargs={}):
    """Plots the orbit data on a specified axis.

    This function takes a DataFrame containing orbit data and creates a plot of the specified columns 
    along the x and y axes. It allows for customization of the plot through additional keyword arguments 
    passed to the `plot` method of Matplotlib's Axes.

    Parameters
    ----------
    data : pandas.DataFrame
        A DataFrame containing the orbit data. The DataFrame must have columns corresponding to the 
        x and y axes specified by `x_axis_column` and `y_axis_column`.

    ax : matplotlib.axes.Axes
        The Matplotlib Axes object where the plot will be drawn. This should be a pre-existing Axes 
        object in which to render the plot.

    x_axis_column : str
        The name of the column in the `data` DataFrame that will be plotted along the x-axis.

    y_axis_column : str
        The name of the column in the `data` DataFrame that will be plotted along the y-axis.

    plotting_kwargs : dict, optional
        A dictionary of additional keyword arguments to customize the appearance of the plot (e.g., line color, 
        linestyle, markers, etc.). These will be passed directly to Matplotlib's `ax.plot()` function.

    Returns
    -------
    None
        This function does not return any value but modifies the `ax` object to include the plotted data.

    Notes
    -----
    - The function assumes that the `data` DataFrame contains numerical data in the specified columns.
    - The `ax` parameter should be an instance of Matplotlib's Axes object. If you do not have an existing Axes object, 
      you can create one using `plt.subplots()` or similar functions in Matplotlib.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>> import pandas as pd
    >>> data = pd.DataFrame({'time': [0, 1, 2, 3], 'x_pos': [0, 1, 2, 3], 'y_pos': [0, 1, 4, 9]})
    >>> fig, ax = plt.subplots()
    >>> orbit_plot(data, ax, x_axis_column='x_pos', y_axis_column='y_pos', plotting_kwargs={'color': 'blue'})
    >>> plt.show()

    In this example, the `x_pos` column will be plotted on the x-axis and the `y_pos` column on the y-axis. 
    The plot will be drawn with blue lines, as specified by the `plotting_kwargs` argument.
    """
    x_axis_data = data[x_axis_column]
    y_axis_data = data[y_axis_column]
    ax.plot(x_axis_data, y_axis_data, **plotting_kwargs)


def spatial_value_plot(data, ax, x_axis_column, science_column, plotting_kwargs={}):
    """Plots scientific values against spatial coordinates on a specified axis.

    This function takes a DataFrame containing spatial and science data and creates a plot of the 
    specified columns along the x and y axes. The plot can be customized using additional keyword 
    arguments passed to the `plot` method of Matplotlib's Axes.

    Parameters
    ----------
    data : pandas.DataFrame
        A DataFrame containing the spatial and science data. The DataFrame must have columns 
        corresponding to the x-axis (spatial data) and y-axis (scientific values) as specified 
        by `x_axis_column` and `science_column`.

    ax : matplotlib.axes.Axes
        The Matplotlib Axes object where the plot will be drawn. This should be a pre-existing Axes 
        object to render the plot.

    x_axis_column : str
        The name of the column in the `data` DataFrame that will be plotted along the x-axis, 
        typically representing spatial coordinates (e.g., time, position).

    science_column : str
        The name of the column in the `data` DataFrame that will be plotted along the y-axis, 
        typically representing scientific values or measurements associated with the spatial data.

    plotting_kwargs : dict, optional
        A dictionary of additional keyword arguments to customize the appearance of the plot (e.g., 
        line color, linestyle, markers, etc.). These will be passed directly to Matplotlib's `ax.plot()` 
        function.

    Returns
    -------
    None
        This function does not return any value but modifies the `ax` object to include the plotted data.

    Notes
    -----
    - The function assumes that the `data` DataFrame contains numerical data in the specified columns.
    - The `ax` parameter should be an instance of Matplotlib's Axes object. If you do not have an existing Axes object, 
      you can create one using `plt.subplots()` or similar functions in Matplotlib.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>> import pandas as pd
    >>> data = pd.DataFrame({'time': [0, 1, 2, 3], 'velocity': [0, 1, 2, 3]})
    >>> fig, ax = plt.subplots()
    >>> spatial_value_plot(data, ax, x_axis_column='time', science_column='velocity', plotting_kwargs={'color': 'red'})
    >>> plt.show()

    In this example, the `time` column will be plotted on the x-axis and the `velocity` column on the y-axis. 
    The plot will be drawn with red lines, as specified by the `plotting_kwargs` argument.
    """
    ax.plot(data[x_axis_column], data[science_column], **plotting_kwargs)