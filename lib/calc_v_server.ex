defmodule CalcVServer do
  @moduledoc """
  Documentation for CalcVServer.
  """
  @doc """
  Hello world.

  ## Examples

      iex> CalcVServer.hello
      :world

  """
  def hello, do: :world
  import IncompressiveKK.Func
  @compile [:native, {:hipe, [:verbose, :o3]}]


  def calcVel kind, velocitys_field, pressure, bc_field, information, name do
    server = :global.whereis_name(name <> "_" <> Atom.to_string(kind))
    send server, {:calc, velocitys_field, pressure, bc_field, information, self}
    receive do
      {simbol, result, ^server} ->
        case simbol do
          :ok ->
            IO.puts "[Info] velocity calculation finished #{inspect DateTime.utc_now}"
            result
          :error ->
            IO.puts "[Error] #{inspect result}"
            raise RuntimeError
        end
    end
  end

  def genCalcServer name, kind do
    pid = spawn(__MODULE__, :calc_server, [kind])
    :global.register_name(name <> "_" <> Atom.to_string(kind), pid)
    IO.puts "[Info] start calc_V_server <#{inspect pid}>"
  end

  def calc_server kind do
    receive do
      {:calc, velocitys_field, pressure, bc_field, information, client} ->
        IO.puts "[Info] velocity calculation <#{inspect client}> #{inspect DateTime.utc_now}"
        try do
          result = deriveVel kind, velocitys_field, pressure, bc_field, information
          IO.puts "[Info] velocity calculation finished <#{inspect client}> #{inspect DateTime.utc_now}"
          send client, {:ok, result, self}
        rescue
          error ->
            send client, {:error, error, self}
        end
    end
    calc_server kind
  end


  def testMain do
    v_field = List.duplicate([0|List.duplicate(1, 400)], 201)
    p_field = List.duplicate([0|List.duplicate(1, 400)], 201)
    bc_field = for j <- 0..200 do
      for i <- 0..400 do
        if 5<i && i<10 && 5<j && j<10 do
          0
        else
          false
        end
      end end
    CalcVServer.genCalcServer "test", :u
    start_time = DateTime.utc_now
    result = CalcVServer.calcVel :u, {v_field, v_field}, p_field, bc_field,
      %{:x_size => 401,
        :y_size => 201,
        :dx => 0.1,
        :dy => 0.1,
        :dt => 0.01,
        :Re => 70}, "test"
    end_time = DateTime.utc_now
    IO.inspect start_time
    IO.inspect end_time
  end

end # CalcVServer
