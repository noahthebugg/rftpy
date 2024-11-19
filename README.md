# G:RFTpy

G:RFTpy is a python library for dealing with
space calculations.

For example the calculation of 
the Semi-Latus Rectum:
$$p = \frac{b^2}{a}$$
or the reduced Vis-Viva equation:
$$v_{r} = \sqrt{ \frac{\mu}{r} }$$

Furthermore you can use a library of
constants for your calculations:
- Earth radius:
$r_e \approx 6371 km$
- Moon graviational constant: 
$\mu_{moon} \approx 4900 \frac{km^3}{s^2}$

## Usage

Below are some examples for usage.

```python
import rftpy as rft

# get Earth specific constants
earth_radius = rft.Earth.radius
earth_mu = rft.Earth.mu

# returns the 1st cosmic velocity
rft.v_cosmic1(earth_mu, earth_radius)

# returns the Semi-Latus Rectum for an ellipsis
rft.semi_latus_rectum(a=100, b=50)

# returns the angular momentum of an object with
# a given orbital radius and velocity
rft.angular_momentum(earth_radius, v=1000)
```

## Contribution

Pull requests are welcome. For major 
changes, please open an issue first to
discuss what you would like to change.

Please make sure to update tests as 
appropriate.

## License

[MIT](https://link)