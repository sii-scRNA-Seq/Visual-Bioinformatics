import { TestBed } from '@angular/core/testing';

import { CurrentDatasetService } from './current-dataset.service';
import { DatasetInfoService } from './dataset-info.service';
import { MockDatasetInfoService } from './mock-dataset-info.service';
import { first, firstValueFrom } from 'rxjs';

describe('CurrentDatasetService', () => {
  let service: CurrentDatasetService;
  let datasetInfoService: DatasetInfoService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      providers: [
        { provide: DatasetInfoService, useClass: MockDatasetInfoService },
      ]
    });
    service = TestBed.inject(CurrentDatasetService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('constructor', () => {
    it('should update the current dataset to the first in datasetInfo if there is at least one', async () => {
      let currentDataset = await firstValueFrom(service.currentDataset.pipe(first()));
      expect(currentDataset.key).toBe('');
      datasetInfoService = TestBed.inject(DatasetInfoService);
      datasetInfoService.setDatasetInfo();
      currentDataset = await firstValueFrom(service.currentDataset.pipe(first()));
      expect(currentDataset.key).toBe('option1');
    });

    it('should keep the current dataset if datasetInfo contains no datasets', async () => {
      let currentDataset = await firstValueFrom(service.currentDataset.pipe(first()));
      expect(currentDataset.key).toBe('');
      datasetInfoService = TestBed.inject(DatasetInfoService);
      currentDataset = await firstValueFrom(service.currentDataset.pipe(first()));
      expect(currentDataset.key).toBe('');
    });
  });

  describe('setCurrentDataset', () => {
    it('should update the current dataset if the new dataset exists', async () => {
      datasetInfoService = TestBed.inject(DatasetInfoService);
      datasetInfoService.setDatasetInfo();
      let currentDataset = await firstValueFrom(service.currentDataset.pipe(first()));
      expect(currentDataset.key).toBe('option1');
      service.setCurrentDataset('option2');
      currentDataset = await firstValueFrom(service.currentDataset.pipe(first()));
      expect(currentDataset.key).toBe('option2');
    });

    it('should keep the current dataset if the new dataset does not exist', async () => {
      datasetInfoService = TestBed.inject(DatasetInfoService);
      datasetInfoService.setDatasetInfo();
      let currentDataset = await firstValueFrom(service.currentDataset.pipe(first()));
      expect(currentDataset.key).toBe('option1');
      service.setCurrentDataset('option3');
      currentDataset = await firstValueFrom(service.currentDataset.pipe(first()));
      expect(currentDataset.key).toBe('option1');
    });
  });
});
