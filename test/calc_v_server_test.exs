defmodule CalcVServerTest do
  use ExUnit.Case
  doctest CalcVServer

  test "calculate next step velocity" do
    v_field = List.duplicate(List.duplicate([0|List.duplicate(1, 7)], 8), 8)
    p_field = List.duplicate(List.duplicate(List.duplicate(1, 8), 8), 8)
    bc_field = for k <- 0..7 do
      for j <- 0..7 do
        for i <- 0..7 do
          if 4==i && 4==j && 4==k do
            0
          else
            nil
          end
        end end end
    CalcVServer.genCalcServer(:u)
    result = CalcVServer.calcVel :u, {v_field, v_field, v_field}, p_field, bc_field,
      %{:x_size => 8,
        :y_size => 8,
        :z_size => 8,
        :dx => 0.1,
        :dy => 0.1,
        :dz => 0.1,
        :dt => 0.01,
        :Re => 70}
    IO.inspect result
  end
end
