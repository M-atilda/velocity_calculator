#file   func.ex
#author mi-na
#date   18/01/14
#brief  calculate next step's velocity field following differential (Kawamura-Kuwahara) scheme
#       except time differential, all differences are calculated as central differences
#         (when the velosity has positive value, artificial viscocities can be regarded as one-side differential)
defmodule IncompressiveKK.Func do
  @compile [:native, {:hipe, [:verbose, :o3]}]
  @type field :: [[float]]

  @spec deriveVel(atom, {field, field}, field, field, map) :: field
  def deriveVel kind, {x_velocity, y_velocity}=velocitys_field, pressure, bc_field,
    %{:x_size => x_size,
      :y_size => y_size}=information do
    velocity = case kind do
                 :u -> x_velocity
                 :v -> y_velocity
               end
    # x_half_size = round((x_size-1) / 2)
    # y_half_size = round((y_size-1) / 2)
    # left_up_task = Task.async(fn -> deriveVelPartially {0..x_half_size, 0..y_half_size}, kind, velocity, velocitys_field, pressure, bc_field, information end)
    # left_down_task = Task.async(fn -> deriveVelPartially {0..x_half_size, (y_half_size+1)..(y_size-1)}, kind, velocity, velocitys_field, pressure, bc_field, information end)
    # right_up_task = Task.async(fn -> deriveVelPartially {(x_half_size+1)..(x_size-1), 0..y_half_size}, kind, velocity, velocitys_field, pressure, bc_field, information end)
    # right_down_task = Task.async(fn -> deriveVelPartially {(x_half_size+1)..(x_size-1), (y_half_size+1)..(y_size-1)}, kind, velocity, velocitys_field, pressure, bc_field, information end)
    # left_up = Task.await left_up_task
    # left_down = Task.await left_down_task
    # right_up = Task.await right_up_task
    # right_down = Task.await right_down_task
    # up_side = Enum.map :lists.zip(left_up, right_up), fn({l, r}) -> List.to_tuple(l ++ r) end
    # down_side = Enum.map :lists.zip(left_down, right_down), fn({l, r}) -> List.to_tuple(l ++ r) end
    # up_side ++ down_side |> List.to_tuple
    deriveVelPartially({0..(x_size-1), 0..(y_size-1)}, kind, velocity, velocitys_field, pressure, bc_field, information)
    |> Enum.map(&(List.to_tuple &1))
    |> List.to_tuple
  end
  @spec deriveVelPartially({any, any}, atom, field, {field, field}, field, field, map) :: field
  defp deriveVelPartially {x_range, y_range}, kind, velocity, velocitys_field, pressure, bc_field,
    %{:dx => dx,
      :dy => dy,
      :x_size => x_size,
      :y_size => y_size,
      :dt => dt,
      :Re => re} do
    dx2 = dx*dx
    dy2 = dy*dy
    dx4 = 4*dx
    dy4 = 4*dy
    for j <- y_range do
      for i <- x_range do
        if !id(bc_field, {i,j}) do
          cond do
            i == (x_size-1) || i == (x_size-2) ->
              2 * id(velocity, {i-1,j}) - id(velocity, {i-2,j})
            j == 0 || j == 1 ->
              2 * id(velocity, {i,j+1}) - id(velocity, {i,j+2})
            j == (y_size-1) || j == (y_size-2) ->
              2 * id(velocity, {i,j-1}) - id(velocity, {i,j-2})
            true ->
              pg_result = calcPreGrad kind, i,j, pressure, dx,dy, x_size,y_size, bc_field
              d_result = calcDiffusion i,j, velocity, dx2,dy2, x_size,y_size, re
              av_result = calcArtiVisc i,j, velocity, velocitys_field, bc_field, dx4,dy4, x_size,y_size
              id(velocity, {i,j}) + dt * (-pg_result + d_result - av_result)
          end
        else
          #TODO: more flexible
          #NOTE: eliminate edges condition
          if 4<i && i<(x_size-1-4) && 4<j && j<(y_size-1-4) do
            cond do
              !id(bc_field, {i-2,j}) && id(bc_field, {i-1,j}) ->
                2 * id(velocity, {i-1,j}) - id(velocity, {i-2,j})
              !id(bc_field, {i+2,j}) && id(bc_field, {i+1,j}) ->
                2 * id(velocity, {i+1,j}) - id(velocity, {i+2,j})
              !id(bc_field, {i,j-2}) && id(bc_field, {i,j-1}) ->
                2 * id(velocity, {i,j-1}) - id(velocity, {i,j-2})
              !id(bc_field, {i,j+2}) && id(bc_field, {i,j+1}) ->
                2 * id(velocity, {i,j+1}) - id(velocity, {i,j+2})
              true ->
                id(bc_field, {i,j})
            end
          else
            id(bc_field, {i,j})
          end
        end
      end
    end
  end


  @spec calcPreGrad(atom, integer, integer, field, float, float, integer, integer, map) :: float
  defp calcPreGrad(:u, i,j, pressure, dx,_dy, x_size,y_size, bc_field) when 0<i and 0<j and i<(x_size-1) and j<(y_size-1) do
    (id(pressure, {i+1,j}) - id(pressure, {i-1,j})) / (2 * dx)
  end
  defp calcPreGrad :u, i,j, pressure, dx,_dy, x_size, _y_size, _bc_field do
    # min_i = max 0, i-1
    # max_i = min (x_size-1), i+1
    # (id(pressure, {max_i,j}) - id(pressure, {min_i,j})) / dx
    0.0
  end
  defp calcPreGrad(:v, i,j, pressure, _dx,dy, x_size,y_size, bc_field) when 0<i and 0<j and i<(x_size-1) and j<(y_size-1) do
    (id(pressure, {i,j+1}) - id(pressure, {i,j-1})) / (2 * dy)
  end
  defp calcPreGrad :v, i,j, pressure, _dx,dy, _x_size, y_size, _bc_field do
    # min_j = max 0, j-1
    # max_j = min (y_size-1), j+1
    # (id(pressure, {i,max_j}) - id(pressure, {i,min_j})) / dy
    0.0
  end

  @spec calcDiffusion(integer, integer, field, float, float, integer, integer, float) :: float
  defp calcDiffusion(i,j, velocity, dx2,dy2, x_size,y_size, re) when 0<i and 0<j and i<(x_size-1) and j<(y_size-1) do
    df2dx2 = calcDiffusionHelper id(velocity, {i+1,j}), id(velocity, {i,j}), id(velocity, {i-1,j}), dx2
    df2dy2 = calcDiffusionHelper id(velocity, {i,j+1}), id(velocity, {i,j}), id(velocity, {i,j-1}), dy2
    (df2dx2 + df2dy2) / re
  end
  defp calcDiffusion i,j, velocity, dx2,dy2, x_size,y_size, re do
    # min_i = max 0, i-1
    # max_i = min (x_size-1), i+1
    # min_j = max 0, j-1
    # max_j = min (y_size-1), j+1
    # df2dx2 = calcDiffusionHelper id(velocity, {max_i,j}), id(velocity, {i,j}), id(velocity, {min_i,j}), dx2
    # df2dy2 = calcDiffusionHelper id(velocity, {i,max_j}), id(velocity, {i,j}), id(velocity, {i,min_j}), dy2
    # (df2dx2 + df2dy2) / re
    0.0
  end
  @spec calcDiffusionHelper(float, float, float, float) :: float
  defp calcDiffusionHelper f_ij1p, f_ij, f_ij1m, delta2 do
    ((f_ij1p - f_ij) - (f_ij - f_ij1m)) / delta2
  end


  @spec calcArtiVisc(integer, integer, field, {field, field}, map, float,float, integer,integer) :: float
  defp calcArtiVisc(i,j, velocity, velocitys_field, bc_field, dx4,dy4, x_size,y_size) when 1<i and 1<j and i<(x_size-2) and j<(y_size-2) do
    #TODO: more flexible (it may not suitable for circle column)
    #NOTE: remove edge's conditions
    udfdx = calcArtiViscX i,j, velocity, velocitys_field, dx4
    vdfdy = calcArtiViscY i,j, velocity, velocitys_field, dy4
    udfdx + vdfdy
  end
  defp calcArtiVisc(_i,_j, _velocity, _velocitys_field, _bc_field, _dx,_dy, _x_size,_y_size), do: 0
  @spec calcArtiViscX(integer, integer, field, {field, field}, float) :: float
  defp calcArtiViscX i,j, f, {u, _v}, dx4 do
    f_ij = id(f, {i,j})
    f_ij1p = id(f, {i+1,j})
    f_ij2p = id(f, {i+2,j})
    f_ij1m = id(f, {i-1,j})
    f_ij2m = id(f, {i-2,j})
    u_ij = id(u, {i,j})
    calcArtiViscHelper u_ij, f_ij, f_ij1p, f_ij2p, f_ij1m, f_ij2m, dx4
  end
  @spec calcArtiViscY(integer, integer, field, {field, field}, float) :: float
  defp calcArtiViscY i,j, f, {_u, v}, dy4 do
    f_ij = id(f, {i,j})
    f_ij1p = id(f, {i,j+1})
    f_ij2p = id(f, {i,j+2})
    f_ij1m = id(f, {i,j-1})
    f_ij2m = id(f, {i,j-2})
    v_ij = id(v, {i,j})
    calcArtiViscHelper v_ij, f_ij, f_ij1p, f_ij2p, f_ij1m, f_ij2m, dy4
  end
  @spec calcArtiViscHelper(float, float, float, float, float, float, float) :: float
  defp calcArtiViscHelper vel_ij, f_ij, f_ij1p, f_ij2p, f_ij1m, f_ij2m, delta4 do
    (vel_ij * ((-f_ij2p + 8*f_ij1p - 8*f_ij1m + f_ij2m) / (3 * delta4)) + abs(vel_ij) * ((f_ij2p - 4*f_ij1p + 6*f_ij - 4*f_ij1m + f_ij2m) / delta4))
  end


  #NOTE: utility functions
  @spec id(field, {integer, integer}) :: any
  def id enumerable, {i, j} do
    elem(elem(enumerable, j), i)
  end

  def getFromKind kind, {x_factor, y_factor} do
    case kind do
      :u ->
        x_factor
      :v ->
        y_factor
    end
  end


end #IncompressiveKK.func
