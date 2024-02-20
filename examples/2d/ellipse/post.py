from numpy import abs, arange, array, exp, log, polyfit, sqrt, zeros
from vtk import vtkUnstructuredGridReader, vtkDataObject, vtkThreshold
from matplotlib.pyplot import subplots, text


def post(mesh_index, mesh_size):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(f'm{mesh_index}/end.vtk')
    reader.Update()

    threshold = vtkThreshold()

    threshold.SetInputConnection(reader.GetOutputPort())
    threshold.SetLowerThreshold(-5.0*mesh_size)
    threshold.SetUpperThreshold( 5.0*mesh_size)
    threshold.SetInputArrayToProcess(
        0, 0, 0,
        vtkDataObject.FIELD_ASSOCIATION_POINTS,
        'level_set')
    threshold.Update()

    grid = threshold.GetOutput()
    grid = reader.GetOutput()
    points = grid.GetPoints()
    level_set = grid.GetPointData().GetArray('level_set')

    difference = zeros(points.GetNumberOfPoints())

    for i in range(points.GetNumberOfPoints()):
        x, y, _ = points.GetPoint(i)
        ls_exact = sqrt(4*x**2 + 9*y**2) - 1.0
        ls_numerical = level_set.GetValue(i)
        difference[i] = abs(ls_exact - ls_numerical)

    return difference.max(), difference.mean()


def fit(x, y):
    m, b = polyfit(x, y, 1)
    return m, b


if __name__ == "__main__":
    exp_range = arange(0, 4)
    mesh_sizes = 0.06 * 2.0**-exp_range

    results = array([post(i, h) for i, h in enumerate(mesh_sizes)])
    fig, axes = subplots(2, 1, sharex=True)


    axes[0].loglog(mesh_sizes, results[:, 0], 's')
    m, b = fit(log(mesh_sizes), log(results[:, 0]))
    axes[0].loglog(mesh_sizes, exp(m*log(mesh_sizes)+b), '-')
    axes[0].text(0.1, 0.5, f"y={m:0.03f}x+{b:0.2f}", transform=axes[0].transAxes)
    axes[1].loglog(mesh_sizes, results[:, 1], 's')
    m, b = fit(log(mesh_sizes), log(results[:, 1]))
    axes[1].loglog(mesh_sizes, exp(m*log(mesh_sizes)+b), '-')
    axes[1].text(0.1, 0.5, f"y={m:0.03f}x+{b:0.2f}", transform=axes[1].transAxes)

    axes[0].set_xlabel('Mesh size')
    axes[0].set_ylabel('Max norm')
    axes[1].set_ylabel('L1 norm')
    for ax in axes:
        ax.grid(True, which='both', ls='--')
    fig.savefig('result.png', bbox_inches='tight')



'''
# Create a vtkPoints object to store X, Y coordinates
    points = vtk.vtkPoints()
    points.SetData(grid.GetPoints())

# Compute function on X, Y coordinates
num_points = points.GetNumberOfPoints()
values = vtk.vtkDoubleArray()
values.SetNumberOfComponents(1)
values.SetName("FunctionValues")
    reader.
'''

