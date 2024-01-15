import { fakeAsync, tick, TestBed } from '@angular/core/testing';
import { first, firstValueFrom, from } from 'rxjs';
import { HttpClientTestingModule, HttpTestingController } from '@angular/common/http/testing';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { Block } from './block.interface';
import { OutputService } from './output.service';

describe('OutputService', () => {
  let service: OutputService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        HttpClientTestingModule,
        MatSnackBarModule,
      ],
    });
    service = TestBed.inject(OutputService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('outputs', () => {
    it('should initially be empty', async () => {
      const outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
    });
  });

  describe('executeBlock', () => {
    it('should add a response to outputs array', fakeAsync(() => {
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
        onRun: () => from(''),
      };
      const mockHttp = TestBed.inject(HttpTestingController);
      service.executeBlock(block).then(async () => {
        const outputs = await firstValueFrom(service.outputs);
        expect(outputs).toEqual([{text: 'Hello World'}]);
      });
      tick();
      const req = mockHttp.expectOne('http://127.0.0.1:5000/loaddata');
      expect(req.request.method).toBe('GET');
      req.flush({text: 'Hello World'});
      mockHttp.verify();  
    }));
  }); 
});
