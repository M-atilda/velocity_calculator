defmodule CalcVServerTest do
  use ExUnit.Case
  doctest CalcVServer
  @tag timeout: 600000

  test "calculate next step velocity" do
    v_field = List.duplicate(List.duplicate([0|List.duplicate(1, 100)], 101), 101)
    p_field = List.duplicate(List.duplicate([0|List.duplicate(1, 100)], 101), 101)
    bc_field = for k <- 0..100 do
      for j <- 0..100 do
        for i <- 0..100 do
          if 24==i || 24==j || 24==k do
            0
          else
            nil
          end
        end end end
    CalcVServer.genCalcServer(:u)
    result = CalcVServer.calcVel :u, {v_field, v_field, v_field}, p_field, bc_field,
      %{:x_size => 101,
        :y_size => 101,
        :z_size => 101,
        :dx => 0.1,
        :dy => 0.1,
        :dz => 0.1,
        :dt => 0.01,
        :Re => 70}
    IO.inspect result
  end
end
