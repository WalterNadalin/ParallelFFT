#### Initial diffusion ###############################
set term png
set output 'plots/diffusivity.png'
min = 0
max = 0.9
set cbrange [min:max]
plot 'data/diffusivity.dat' matrix with image

#### Initial concentration ###########################
set term png
set output 'plots/initial_concentration.png'
min = 0
max = 0.9
set cbrange [min:max]
plot 'data/initial_concentration.dat' matrix with image

#### Concentration evolving in time ###################
file_exists(file) = system("[ -f '".file."' ] && echo '1' || echo '0'") + 0
set term gif animate
unset key
set title 'Concentration in time'
set output 'plots/animation.gif'
frames = 17
min = 0
max = 0.5
set cbrange [min: max]
i = 1

while (file_exists('data/concentration_'.i.'.dat')) {
  plot 'data/concentration_'.i.'.dat' matrix  with image
  i = i + 1
}

