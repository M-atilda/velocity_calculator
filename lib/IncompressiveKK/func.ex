#file   func.ex
#author mi-na
#date   18/01/14
#brief  calculate next step's velocity field following differential (Kawamura-Kuwahara) scheme
#       except time differential, all differences are calculated as central differences
#         (when the velosity has positive value, artificial viscocities can be regarded as one-side differential)
defmodule IncompressiveKK.Func do
  import IncompressiveKK.Util


  def deriveVel kind, {x_velocity, y_velocity}=velocitys_field, pressure, bc_field,
    %{:x_size => x_size,
      :y_size => y_size,
      :dt => dt}=information do
    velocity = case kind do
                 :u -> x_velocity
                 :v -> y_velocity
               end
    pressure_gradient = Task.async(__MODULE__, :calcPreGrad, [kind, pressure, information])
    diffusion = Task.async(__MODULE__, :calcDiffusion, [velocity, information])
    artificial_viscocity = Task.async(__MODULE__, :calcArtiVisc, [velocity, velocitys_field, information])
    pg_result = Task.await(pressure_gradient)
    d_result = Task.await(diffusion)
    av_result = Task.await(artificial_viscocity)
    for j <- 0..(y_size-1) do
      for i <- 0..(x_size-1) do
        if id(bc_field, {i,j}) == nil do
          id(velocity, {i,j}) + dt * (-id(pg_result, {i,j}) + id(d_result, {i,j}) - id(av_result, {i,j}))
        else
          id(bc_field, {i,j})
        end
      end
    end
  end


  def calcPreGrad :u, pressure, %{:dx => dx,
                                  :x_size => x_size, :y_size => y_size} do
    for j <- 0..(y_size-1) do
      for i <- 0..(x_size-1) do
        if 0<i && 0<j && i<(x_size-1) && j<(y_size-1) do
          (id(pressure, {i+1,j}) - id(pressure, {i-1,j})) / (2 * dx)
        else
          min_i = max 0, i-1
          max_i = min (x_size-1), i+1
          (id(pressure, {max_i,j}) - id(pressure, {min_i,j})) / dx
        end
      end end
  end
  def calcPreGrad :v, pressure, %{:dy => dy, :x_size => x_size, :y_size => y_size} do
    for j <- 0..(y_size-1) do
      for i <- 0..(x_size-1) do
        if 0<i && 0<j && i<(x_size-1) && j<(y_size-1) do
          (id(pressure, {i,j+1}) - id(pressure, {i,j-1})) / (2 * dy)
        else
          min_j = max 0, j-1
          max_j = min (y_size-1), j+1
          (id(pressure, {i,max_j}) - id(pressure, {i,min_j})) / dy
        end
      end end
  end


  def calcDiffusion velocity, %{:dx => dx,
                                :dy => dy,
                                :x_size => x_size,
                                :y_size => y_size,
                                :Re => re} do
    for j <- 0..(y_size-1) do
      for i <- 0..(x_size-1) do
        if 0<i && 0<j && i<(x_size-1) && j<(y_size-1) do
          df2dx2 = calcDiffusionHelper id(velocity, {i+1,j}), id(velocity, {i,j}), id(velocity, {i-1,j}), dx
          df2dy2 = calcDiffusionHelper id(velocity, {i,j+1}), id(velocity, {i,j}), id(velocity, {i,j-1}), dy
          (df2dx2 + df2dy2) / re
        else
          min_i = max 0, i-1
          max_i = min (x_size-1), i+1
          min_j = max 0, j-1
          max_j = min (y_size-1), j+1
          df2dx2 = calcDiffusionHelper id(velocity, {max_i,j}), id(velocity, {i,j}), id(velocity, {min_i,j}), dx
          df2dy2 = calcDiffusionHelper id(velocity, {i,max_j}), id(velocity, {i,j}), id(velocity, {i,min_j}), dy
          (df2dx2 + df2dy2) / re
        end
      end
    end
  end
  defp calcDiffusionHelper f_ij1p, f_ij, f_ij1m, delta do
    ((f_ij1p - f_ij) - (f_ij - f_ij1m)) / (delta * delta)
  end
  

  def calcArtiVisc velocity, velocitys_field, %{:dx => dx,
                                                :dy => dy,
                                                :x_size => x_size,
                                                :y_size => y_size} do
    for j <- 0..(y_size-1) do
      for i <- 0..(x_size-1) do
        if 1<i && 1<j && i<(x_size-2) && j<(y_size-2) do
          udfdx = calcArtiViscX {i,j}, velocity, velocitys_field, dx
          vdfdy = calcArtiViscY {i,j}, velocity, velocitys_field, dy
          udfdx + vdfdy
        else
          0
        end
      end
    end
  end
  defp calcArtiViscX {i,j}, f, {u, _v}, dx do
    f_ij = id(f, {i,j})
    f_ij1p = id(f, {i+1,j})
    f_ij2p = id(f, {i+2,j})
    f_ij1m = id(f, {i-1,j})
    f_ij2m = id(f, {i-2,j})
    u_ij = id(u, {i,j})
    calcArtiViscHelper u_ij, f_ij, f_ij1p, f_ij2p, f_ij1m, f_ij2m, dx
  end
  defp calcArtiViscY {i,j}, f, {_u, v}, dy do
    f_ij = id(f, {i,j})
    f_ij1p = id(f, {i,j+1})
    f_ij2p = id(f, {i,j+2})
    f_ij1m = id(f, {i,j-1})
    f_ij2m = id(f, {i,j-2})
    v_ij = id(v, {i,j})
    calcArtiViscHelper v_ij, f_ij, f_ij1p, f_ij2p, f_ij1m, f_ij2m, dy
  end
  defp calcArtiViscHelper vel_ij, f_ij, f_ij1p, f_ij2p, f_ij1m, f_ij2m, delta do
    vel_ij * ((-f_ij2p + 8*f_ij1p - 8*f_ij1m + f_ij2m) / (12 * delta)) + abs(vel_ij) * ((f_ij2p - 4*f_ij1p + 6*f_ij - 4*f_ij1m + f_ij2m) / (4 * delta))
  end


end #IncompressiveKK.func
