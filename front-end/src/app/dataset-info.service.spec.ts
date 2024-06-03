import { TestBed } from '@angular/core/testing';

import { DatasetInfoService } from './dataset-info.service';

describe('DatasetInfoService', () => {
  let service: DatasetInfoService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: []
    });
    service = TestBed.inject(DatasetInfoService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
