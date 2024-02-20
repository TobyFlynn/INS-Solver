from numpy import abs, arange, array, exp, log, polyfit, sqrt, zeros
from vtk import vtkUnstructuredGridReader, vtkThreshold, vtkDataObject
from matplotlib.pyplot import subplots, text


def post(mesh_index, mesh_size):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(f'm{mesh_index}/end.vtk')
    reader.Update()

    threshold = vtkThreshold()
    threshold.SetInputConnection(reader.GetOutputPort())
    threshold.SetThresholdFunction(0)
    threshold.SetLowerThreshold(-5.0*mesh_size)
    threshold.SetUpperThreshold( 5.0*mesh_size)
    threshold.SetInputArrayToProcess(
        0, 0, 0,
        vtkDataObject.FIELD_ASSOCIATION_POINTS,
        'level_set')
    threshold.Update()

    grid = threshold.GetOutput()
    points = grid.GetPoints()
    level_set = grid.GetPointData().GetArray('level_set')

    difference = zeros(points.GetNumberOfPoints())

    for i in range(points.GetNumberOfPoints()):
        x, y, _ = points.GetPoint(i)
        ls_exact = sqrt(x**2 + y**2) - 1.0
        ls_numerical = level_set.GetValue(i)
        difference[i] = abs(ls_exact - ls_numerical)

    return difference.max(), difference.mean()


def fit(x, y):
    m, b = polyfit(x, y, 1)
    return m, b


if __name__ == "__main__":
    exp_range = arange(0, 4)
    mesh_sizes = 0.06 * 2.0**-exp_range

    # results = array([post(i, h) for i, h in enumerate(mesh_sizes)])
    results = array([[
        4.3451e-06, 2.92192e-07, 1.79148e-08, 1.08719e-09],
        [1.2085e-05, 1.54777e-06, 6.05132e-08, 3.53746e-09]]).T

    fig, axes = subplots(2, 1, sharex=True)


    for i in range(2):
        axes[i].loglog(mesh_sizes, results[:, i], 's')
        m, b = fit(log(mesh_sizes), log(results[:, i]))
        axes[i].loglog(mesh_sizes, exp(m*log(mesh_sizes)+b), '-')
        axes[i].text(0.1, 0.5, f"y={m:0.02f}x+{b:0.2f}", transform=axes[i].transAxes)

    axes[0].set_xlabel('Mesh size')
    axes[0].set_ylabel('Max norm')
    axes[1].set_ylabel('L1 norm')
    for ax in axes:
        ax.grid(True, which='both', ls='--')
    fig.savefig('result.png', bbox_inches='tight')
