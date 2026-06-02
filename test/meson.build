tests = [
	'TESTF.f',
	'TESTP.f',
	'TESTS.f',
]

foreach t : tests
  test(
    t,
    executable(
      'test-@0@'.format(t.replace('.f', '')),
      sources: [t],
      link_with: hompack_lib,
    ),
    suite: 'unit',
    workdir: meson.project_source_root(),
  )
endforeach
