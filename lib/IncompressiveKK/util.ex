defmodule IncompressiveKK.Util do

  def id enumerable, {i, j, k} do
    Enum.at(Enum.at(Enum.at(enumerable, k), j), i)
  end

  def getFromKind kind, {x_factor, y_factor, z_factor} do
    case kind do
      :u ->
        x_factor
      :v ->
        y_factor
      :w ->
        z_factor
    end
  end

end # IncompressiveKK.Util
