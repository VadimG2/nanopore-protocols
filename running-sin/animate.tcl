# 1. Загрузка структуры и траектории
mol new scaled_sin_ions.pdb
mol addfile sin_20V_restart.dcd waitfor all

# 2. Настройка вида
rotate x by 90
display projection orthographic
display background white
axes location off

# Важное обновление вида – центрирование молекул
display resetview

# 3. Масштабирование и обновление отображения
scale by 1.5
display update

# 4. Подготовка анимации
animate goto start
set numframes [molinfo top get numframes]

# Создаем каталог для сохранения кадров, если он не существует
if {![file isdirectory "movie/first"]} {
    file mkdir "movie/first"
}

# 5. Рендеринг кадров
for {set i 0} {$i < $numframes} {incr i} {
    animate goto $i
    display update
    render Snapshot -format ppm [format "movie/first/frame%04d.ppm" $i]
}

puts "Сохранение кадров завершено!"
