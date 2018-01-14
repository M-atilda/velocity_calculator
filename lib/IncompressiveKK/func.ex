#file   func.ex
#author mi-na
#date   18/01/14
#brief  calculate next step's velocity field following differential (Kawamura-Kuwahara) scheme
#       except time differential, all differences are calculated as central differences
#         (when the velosity has positive value, artificial viscocities can be regarded as one-side differential)
defmodule IncompressiveKK.Func do
  import IncompressiveKK.Util


  def deriveVel kind, {x_velocity, y_velocity, z_velocity}=velocitys_field, pressure, bc_field,
    %{:x_size => x_size,
      :y_size => y_size,
      :z_size => z_size,
      :dt => dt}=information do
    velocity = case kind do
                 :u -> x_velocity
                 :v -> y_velocity
                 :w -> z_velocity
               end
    pressure_gradient = Task.async(calcPreGrad kind, pressure, information)
    diffusion = Task.async(calcDiffusion velocity, information)
    artificial_viscocity = Task.async(calcArtiVisc velocity, velocitys_field, information)
    pg_result = Task.await(pressure_gradient)
    d_result = Task.await(diffusion)
    av_result = Task.await(artificial_viscocity)
    for k <- 0..z_size do
      for j <- 0..y_size do
        for i <- 0..x_size do
          if id(bc_field, {i,j,k}) == nil do
            id(velocity, {i,j,k}) +
            dt * ((-id(pg_result, {i,j,k})) +
              id(d_result, {i,j,k}) -
              id(av_result, {i,j,k}))
          else
            id(bc_field, {i,j,k})
          end
        end
      end
    end
  end


  def calcPreGrad :u, pressure, %{:dx => dx,
                                  :x_size => x_size, :y_size => y_size, :z_size => z_size} do
    for k <- 0..z_size do
      for j <- 0..y_size do
        for i <- 0..x_size do
            if 0<i && 0<j && 0<k && i<x_size && j<y_size && k<z_size do
              (id(pressure, {i+1,j,k}) - id(pressure, {i-1,j,k})) / (2 * dx)
            else
              #NOTE: should these differentials are calulated using one-side differential?
              0
            end
        end end end
  end
  def calcPreGrad :v, pressure, %{:dy => dy, :x_size => x_size, :y_size => y_size, :z_size => z_size} do
    for k <- 0..z_size do
      for j <- 0..y_size-1 do
        for i <- 0..x_size do
          if 0<i && 0<j && 0<k && i<x_size && j<y_size && k<z_size do
            (id(pressure, {i,j+1,k}) - id(pressure, {i,j-1,k})) / (2 * dy)
          else
            0
          end
        end end end
  end
  def calcPreGrad :w, pressure, %{:dz => dz, :x_size => x_size, :y_size => y_size, :z_size => z_size} do
    for k <- 0..z_size do
      for j <- 0..y_size do
        for i <- 0..x_size do
          if 0<i && 0<j && 0<k && i<x_size && j<y_size && k<z_size do
            (id(pressure, {i,j,k+1}) - id(pressure, {i-1,j,k-1})) / (2 * dz)
          else
            0
          end
        end end end
  end


  def calcDiffusion velocity, %{:dx => dx,
                                :dy => dy,
                                :dz => dz,
                                :x_size => x_size,
                                :y_size => y_size,
                                :z_size => z_size,
                                :Re => re} do
    for k <- 0..z_size do
      for j <- 0..y_size do
        for i <- 0..x_size do
          if 1<i && 1<j && 1<k && i<(x_size-1) && j<(y_size-1) && k<(z_size-1) do
            df2dx2 = calcDiffusionHelper id(i+1,j,k), id(i,j,k), id(i-1,j,k), dx
            df2dy2 = calcDiffusionHelper id(i,j+1,k), id(i,j,k), id(i,j-1,k), dy
            df2dz2 = calcDiffusionHelper id(i,j,k+1), id(i,j,k), id(i,j,k-1), dz
            (df2dx2 + df2dy2 + df2dz2) / re
          else
            0
          end
        end
      end
    end
  end
  defp calcDiffusionHelper f_ijk1p, f_ijk, f_ijk1m, delta do
    ((f_ijk1p - f_ijk) - (f_ijk - f_ijk1m)) / (delta * delta)
  end
  

  def calcArtiVisc velocity, velocitys_field, %{:dx => dx,
                                                :dy => dy,
                                                :dz => dz,
                                                :x_size => x_size,
                                                :y_size => y_size,
                                                :z_size => z_size} do
    for k <- 0..z_size do
      for j <- 0..y_size do
        for i <- 0..x_size do
          if 1<i && 1<j && 1<k && i<(x_size-1) && j<(y_size-1) && k<(z_size-1) do
            udfdx = calcArtiViscX {i,j,k}, velocity, velocitys_field, dx
            vdfdy = calcArtiViscY {i,j,k}, velocity, velocitys_field, dy
            wdfdz = calcArtiViscZ {i,j,k}, velocity, velocitys_field, dz
            udfdx + vdfdy + wdfdz
          else
            0
          end
        end
      end
    end
  end
  defp calcArtiViscX {i,j,k}, f, {u, _v, _w}, dx do
    f_ijk = id(f, {i, j, k})
    f_ijk1p = id(f, {i+1,j,k})
    f_ijk2p = id(f, {i+2,j,k})
    f_ijk1m = id(f, {i-1,j,k})
    f_ijk2m = id(f, {i-2,j,k})
    u_ijk = id(u, {i, j, k)
    calcArtiViscHelper u_ijk, f_ijk, f_ijk1p, f_ijk2p, f_ijk1m, f_ijk2m, dx
  end
  defp calcArtiViscY pos, f, {_u, v, _w}, dy do
    f_ijk = id(f, {i,j,k})
    f_ijk1p = id(f, {i,j+1,k})
    f_ijk2p = id(f, {i,j+2,k})
    f_ijk1m = id(f, {i,j-1,k})
    f_ijk2m = id(f, {i,j-2,k})
    v_ijk = id(v, {i,j,k)
    calcArtiViscHelper v_ijk, f_ijk, f_ijk1p, f_ijk2p, f_ijk1m, f_ijk2m, dy
  end
  defp calcArtiViscZ pos, f, {_u, _v, w}, dz do
    f_ijk = id(f, {i,j,k})
    f_ijk1p = id(f, {i,j,k+1})
    f_ijk2p = id(f, {i,j,k+2})
    f_ijk1m = id(f, {i,j,k-1})
    f_ijk2m = id(f, {i,j,k-2})
    w_ijk = id(w, {i,j,k)
    calcArtiViscHelper w_ijk, f_ijk, f_ijk1p, f_ijk2p, f_ijk1m, f_ijk2m, dz
  end
  defp calcArtiViscHelper vel_ijk, f_ijk, f_ijk1p, f_ijk2p, f_ijk1m, f_ijk2m, delta do
    vel_ijk *
    ((-f_ijk2p + 8*f_ijk1p - 8*f_ijk1m + f_ijk2m) / (12 * delta)) +
    abs(vel_ijk) *
    ((f_ijk2p - 4*f_ijk1p + 6*f_ijk - 4*f_ijk1m + f_ijk2m) / (4 * delta))    
  end



end #IncompressiveKK.func
