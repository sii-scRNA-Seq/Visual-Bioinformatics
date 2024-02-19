import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { fakeAsync, tick, TestBed } from '@angular/core/testing';
import { first, firstValueFrom, from } from 'rxjs';
import { HttpClientTestingModule, HttpTestingController } from '@angular/common/http/testing';
import { MatSnackBar, MatSnackBarModule } from '@angular/material/snack-bar';

import { Block } from './block.interface';
import { MockUserIdService } from './mock-user-id.service';
import { OutputService } from './output.service';
import { UserIdService } from './user-id.service';

describe('OutputService', () => {
  let service: OutputService;
  let snackBar: MatSnackBar;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        BrowserAnimationsModule,
        HttpClientTestingModule,
        MatSnackBarModule,
      ], 
      providers: [
        { provide: UserIdService, useClass: MockUserIdService }
      ],
    });
    service = TestBed.inject(OutputService);
    snackBar = TestBed.inject(MatSnackBar);
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
    it('should add a response to outputs array when it receives a valid response', fakeAsync(() => {
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
        onRun: () => from(''),
      };
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const mockHttp = TestBed.inject(HttpTestingController);
      service.executeBlock(block).then(async () => {
        const outputs = await firstValueFrom(service.outputs);
        expect(outputs).toEqual([{text: 'Hello World'}]);
      });
      tick();
      const req = mockHttp.expectOne('http://127.0.0.1:5000/loaddata?user_id=mock_user_id');
      expect(req.request.method).toBe('GET');
      req.flush({text: 'Hello World'});
      mockHttp.verify();
    }));

    it('should open snack bar when userId is null', fakeAsync(() => {
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
        onRun: () => from(''),
      };
      const spy = spyOn(snackBar, 'open');
      service.executeBlock(block).then(async () => {
        expect(spy).toHaveBeenCalledOnceWith('No User ID, please refresh the page and try again', 'Close', { duration: 5000 });
      });
    }));

    it('should open snack bar when response is an error', fakeAsync(() => {
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
        onRun: () => from(''),
      };
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const mockHttp = TestBed.inject(HttpTestingController);
      const spy = spyOn(snackBar, 'open');
      service.executeBlock(block).then(async () => {
        expect(spy).toHaveBeenCalledOnceWith('Not a valid order, please check the blocks and try again', 'Close', { duration: 5000 });
      });
      tick();
      const req = mockHttp.expectOne('http://127.0.0.1:5000/loaddata?user_id=mock_user_id');
      expect(req.request.method).toBe('GET');
      req.flush('', { status: 406, statusText: 'Bad Request'});
      mockHttp.verify(); 
    }));
  }); 
});
