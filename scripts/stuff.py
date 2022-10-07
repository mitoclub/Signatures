from yellowbrick.datasets import load_concrete
from yellowbrick.features import PCA
from yellowbrick.exceptions import YellowbrickValueError


class MyPCA(PCA):
    def _draw_projection_features(self, Xp, y):
        x_vector = self.pca_components_[0]
        y_vector = self.pca_components_[1]
        max_x = max(Xp[:, 0])
        max_y = max(Xp[:, 1])
        if self.projection == 2:
            for i in range(self.pca_components_.shape[1]):
                self.ax.arrow(
                    x=0,
                    y=0,
                    dx=x_vector[i] * max_x,
                    dy=y_vector[i] * max_y,
                    color="r",
                    head_width=0.05,
                    width=0.005,
                )
                self.ax.text(
                    x_vector[i] * max_x * 1.05,
                    y_vector[i] * max_y * 1.05,
                    self.features_[i],
                    color="r",
                )
        elif self.projection == 3:
            z_vector = self.pca_components_[2]
            max_z = max(Xp[:, 1])
            for i in range(self.pca_components_.shape[1]):
                self.ax.plot(
                    [0, x_vector[i] * max_x],
                    [0, y_vector[i] * max_y],
                    [0, z_vector[i] * max_z],
                    color="r",
                )
                self.ax.text(
                    x_vector[i] * max_x * 1.05,
                    y_vector[i] * max_y * 1.05,
                    z_vector[i] * max_z * 1.05,
                    self.features_[i],
                    color="r",
                )
        else:
            raise YellowbrickValueError("Projection dimensions must be either 2 or 3")

        return self.ax



# Load the concrete dataset
X, y = load_concrete()
print(X.shape, y.shape)

visualizer = MyPCA(scale=True, proj_features=True)
visualizer.fit_transform(X, y)
visualizer.show()
