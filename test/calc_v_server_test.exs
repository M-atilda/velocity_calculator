defmodule CalcVServerTest do
  use ExUnit.Case
  doctest CalcVServer
  @tag timeout: 600000

  test "calculate next step velocity" do
    line = List.to_tuple [0.0|List.duplicate(1.0, 400)]
    v_field = Tuple.duplicate(line, 201)
    p_field = Tuple.duplicate(line, 201)
    bc_field = for j <- 0..200 do
      for i <- 0..400 do
        if 5<i && i<10 && 5<j && j<10 do
          0
        else
          false
        end
      end
      |> List.to_tuple
    end
    |> List.to_tuple
    CalcVServer.genCalcServer(:u)
    start_time = DateTime.utc_now
    result = CalcVServer.calcVel :u, {v_field, v_field}, p_field, bc_field,
      %{:x_size => 401,
        :y_size => 201,
        :dx => 0.1,
        :dy => 0.1,
        :dt => 0.01,
        :Re => 70}
    end_time = DateTime.utc_now
    IO.inspect start_time
    IO.inspect end_time
    IO.inspect result
  end
end
